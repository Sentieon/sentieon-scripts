#!/bin/sh
#
# The script performs germline variant calling from PacBio HiFi reads with Sentieon DNAscope

#### Argument Parsing ####
die()
{
    _ret="${2:-1}"
    test "${_PRINT_HELP:-no}" = yes && print_help >&2
    echo "$1" >&2
    exit "${_ret}"
}


# THE DEFAULTS INITIALIZATION - POSITIONALS
# THE DEFAULTS INITIALIZATION - OPTIONALS
_arg_reference_fasta=
_arg_sample_input=
_arg_model_bundle=
_arg_dbsnp=
_arg_regions=
_arg_threads=
_arg_gvcf="off"


print_help()
{
    printf '%s\n' "Call small variants from PacBio HiFi reads with Sentieon DNAscope"
    printf 'Usage: %s -r <reference> -i <sample-input> -m <model-file> [-d <dbSNP>] [-b <regions>] [-t <threads>] [-g] [-h] [--] <output-vcf>\n' "$0"
    printf '\t%s\n' "<output-vcf>: The output VCF file."
    printf '\t%s\n' "-r: The reference FASTA file. (required)"
    printf '\t%s\n' "-i: The aligned BAM or CRAM file. (required)"
    printf '\t%s\n' "-m: The model bundle file. (required)"
    printf '\t%s\n' "-d: dbSNP VCF file. Supplying this file will annotate variants with their dbSNP refSNP ID numbers. (no default)"
    printf '\t%s\n' "-b: Region BED file. Supplying this file will limit variant calling to the intervals inside the BED file. (no default)"
    printf '\t%s\n' "-t: Number of threads/processes to use. (number of available cores)"
    printf '\t%s\n' "-g: Generate a gVCF output file along with the VCF. (default generates only the VCF)"
    printf '\t%s\n' "-h: Print this help message"
}


parse_commandline()
{
    while getopts 'r:i:m:d:b:t:gh' _key
    do
        case "$_key" in
            r)
                test "x$OPTARG" = x && die "Missing value for the argument '-$_key'." 1
                _arg_reference_fasta="$OPTARG"
                ;;
            i)
                test "x$OPTARG" = x && die "Missing value for the argument '-$_key'." 1
                _arg_sample_input="$OPTARG"
                ;;
            m)
                test "x$OPTARG" = x && die "Missing value for the argument '-$_key'." 1
                _arg_model_bundle="$OPTARG"
                ;;
            d)
                test "x$OPTARG" = x && die "Missing value for the optional argument '-$_key'." 1
                _arg_dbsnp="$OPTARG"
                ;;
            b)
                test "x$OPTARG" = x && die "Missing value for the optional argument '-$_key'." 1
                _arg_regions="$OPTARG"
                ;;
            t)
                test "x$OPTARG" = x && die "Missing value for the optional argument '-$_key'." 1
                _arg_threads="$OPTARG"
                ;;
            g)
                _arg_gvcf="on"
                ;;
            h)
                print_help
                exit 0
                ;;
            *)
                _PRINT_HELP=yes die "FATAL ERROR: Got an unexpected option '-${_key}'" 1
                ;;
        esac
    done
}


handle_passed_args_count()
{
    _required_args_string="'output-vcf'"
    test "${_positionals_count}" -ge 1 || _PRINT_HELP=yes die "ERROR: Not enough positional arguments - we require exactly 1 (namely: $_required_args_string), but got only ${_positionals_count}." 1
    test "${_positionals_count}" -le 1 || _PRINT_HELP=yes die "ERROR: There were spurious positional arguments --- we expect exactly 1 (namely: $_required_args_string), but got ${_positionals_count} (the last one was: '${_last_positional}')." 1
}


assign_positional_args()
{
    _shift_for=$1
    _positional_names="_arg_output_vcf "

    shift "$_shift_for"
    for _positional_name in ${_positional_names}
    do
        test $# -gt 0 || break
        eval "$_positional_name=\${1}" || die "Error during argument parsing, possibly an Argbash bug." 1
        shift
    done
}

parse_commandline "$@"
_positionals_count=$(($# - OPTIND + 1)); _last_positional=$(eval "printf '%s' \"\$$#\""); handle_passed_args_count
assign_positional_args "$OPTIND" "$@"
_positionals_count=$(($# - OPTIND + 1)); _last_positional=$(eval "printf '%s' \"\$$#\""); handle_passed_args_count

case "$_arg_output_vcf" in
    *.vcf.gz)
        ;;
    *)
        die "The output VCF file is expected to have a '.vcf.gz' suffix" 1
        ;;
esac

if [ -z "$_arg_threads" ]; then
    _arg_threads=$(nproc)
fi

if [ -z "$_arg_reference_fasta" ]; then
	_PRINT_HELP=yes die "Error: the '-r' argument is required" 1
fi

if [ -z "$_arg_sample_input" ]; then
	_PRINT_HELP=yes die "Error: the '-i' argument is required" 1
fi

if [ -z "$_arg_model_bundle" ]; then
	_PRINT_HELP=yes die "Error: the '-m' argument is required" 1
fi

#### Check the PATH for the required executable files ####
for exec_file in bcftools bedtools sentieon; do
    if ! which $exec_file 1>/dev/null 2>/dev/null; then
        echo "Error: no '$exec_file'  executable found in the PATH"
        exit 2
    fi
done

#### Check executables have the required versions ####
VERSION_PATTERN='^v\{0,1\}\([[:digit:]]\{1,\}\).\{0,1\}\([[:digit:]]\{1,\}\)\{0,1\}.\{0,1\}\([[:digit:]]\{1,\}\)\{0,1\}$'

cmp() {
    num_a="$1"
    num_b="$2"
    if [ "$num_a" -lt "$num_b" ]; then
        return 1
    elif [ "$num_a" -gt "$num_b" ]; then
        return 2
    fi
    return 0
}

cmp_vers() {
    if [ "$1" = "$2" ];
    then
        return 0
    fi

    a_major=$(echo $1 | sed -e 's/'"$VERSION_PATTERN"'/\1/')
    a_minor=$(echo $1 | sed -e 's/'"$VERSION_PATTERN"'/\2/')
    [ -n "$a_minor" ] || a_minor=0
    a_patch=$(echo $1 | sed -e 's/'"$VERSION_PATTERN"'/\3/')
    [ -n "$a_patch" ] || a_patch=0
    b_major=$(echo $2 | sed -e 's/'"$VERSION_PATTERN"'/\1/')
    b_minor=$(echo $2 | sed -e 's/'"$VERSION_PATTERN"'/\2/')
    [ -n "$b_minor" ] || b_minor=0
    b_patch=$(echo $2 | sed -e 's/'"$VERSION_PATTERN"'/\3/')
    [ -n "$b_patch" ] || b_patch=0

    if cmp "$a_major" "$b_major"; then
        if cmp "$a_minor" "$b_minor"; then
            cmp "$a_patch" "$b_patch"
            return $?
        else
            cmp "$a_minor" "$b_minor"
            return $?
        fi
    else
        cmp "$a_major" "$b_major"
        return $?
    fi
}

SENTIEON_MIN_VERSION=202308
BCFTOOLS_MIN_VERSION=1.10
sentieon_version_string=$(sentieon driver --version)
sentieon_version_string="${sentieon_version_string##sentieon-genomics-}"
bcftools_version_string=$(bcftools --version | head -n 1)
bcftools_version_string="${bcftools_version_string##bcftools }"

if echo "$sentieon_version_string" | grep -q "$VERSION_PATTERN"; then
    cmp_vers "$SENTIEON_MIN_VERSION" "$sentieon_version_string"
    if [ "$?" -eq "2" ]; then
        echo "Error: the pipeline requires sentieon version '$SENTIEON_MIN_VERSION' or later but sentieon '$sentieon_version_string' was found in the PATH"
        exit 2
    fi
else
    echo "Warning: unable to check sentieon version string"
fi

if echo "$bcftools_version_string" | grep -q "$VERSION_PATTERN"; then
    cmp_vers "$BCFTOOLS_MIN_VERSION" "$bcftools_version_string"
    if [ "$?" -eq "2" ]; then
        echo "Error: the pipeline requires bcftools version '$BCFTOOLS_MIN_VERSION' or later but bcftools '$bcftools_version_string' was found in the PATH"
        exit 2
    fi
else
    echo "Warning: unable to check bcftools version string"
fi

#### Check that all input files can be found and are not empty
if [ ! -f "$_arg_reference_fasta" ]; then
    echo "Error: cannot find the input reference fasta file, '$_arg_reference_fasta'"
    exit 2
elif [ ! -s "$_arg_reference_fasta" ]; then
    echo "Error: the input reference fasta file is empty, '$_arg_reference_fasta'"
    exit 2
fi

if [ ! -f "$_arg_sample_input" ]; then
    echo "Error: cannot find the input BAM file, '$_arg_sample_input'"
    exit 2
elif [ ! -s "$_arg_sample_input" ]; then
    echo "Error: the input BAM file is empty, '$_arg_sample_input'"
    exit 2
fi

if [ ! -f "$_arg_model_bundle" ]; then
    echo "Error: cannot find the input model bundle file, '$_arg_model_bundle'"
    exit 2
elif [ ! -s "$_arg_model_bundle" ]; then
    echo "Error: the input model bundle file is empty, '$_arg_model_bundle'"
    exit 2
fi

if [ -n "$_arg_dbsnp" ]; then
    if [ ! -f "$_arg_dbsnp" ]; then
        echo "Error: cannot find the optional dbSNP VCF file, '$_arg_dbsnp'"
        exit 2
    elif [ ! -s "$_arg_dbsnp" ]; then
        echo "Error: the optional dbSNP VCF file is empty, '$_arg_dbsnp'"
        exit 2
    fi
fi

if [ -n "$_arg_regions" ]; then
    if [ ! -f "$_arg_regions" ]; then
        echo "Error: cannot find the optional regions BED file, '$_arg_regions'"
        exit 2
    elif [ ! -s "$_arg_regions" ]; then
        echo "Error: the optional regions BED file is empty, '$_arg_regions'"
        exit 2
    fi
fi

#### Temporary directory ####
tmp_base="/tmp"
if [ -n "$SENTIEON_TMPDIR" ]; then
    tmp_base="$SENTIEON_TMPDIR"
elif [ -n "$TMPDIR" ]; then
    tmp_base="$TMPDIR"
fi
tmp_base="$tmp_base/$$"

tmp_dir="$tmp_base"
while [ -e "$tmp_dir" ]; do
    rand_num=$(echo "" | awk 'BEGIN {srand()} {print int(rand() * 65535)}')
    tmp_dir="${tmp_base}"_"$rand_num"
done
mkdir -p "$tmp_dir"


#### Output files ####
OUTPUT_BASENAME=$(basename "${_arg_output_vcf%%.vcf.gz}")
TMP_BASE="$tmp_dir"/"$OUTPUT_BASENAME"

#### - Diploid
DIPLOID_TMP_OUT="${TMP_BASE}_diploid_tmp.vcf.gz"
DIPLOID_OUT="${TMP_BASE}_diploid.vcf.gz"
DIPLOID_GVCF="${TMP_BASE}_diploid.g.vcf.gz"

#### - Phasing
PHASED_VCF="${TMP_BASE}_diploid_phased.vcf.gz"
PHASED_BED="${TMP_BASE}_diploid_phased.bed"
PHASED_EXT="${TMP_BASE}_diploid_phased.ext.vcf.gz"
REF_BED="${TMP_BASE}_reference.bed"
UNPHASED_BED="${TMP_BASE}_diploid_unphased.bed"
PHASED_UNPHASED_VCF="${TMP_BASE}_diploid_phased_unphased.vcf.gz"
REPEAT_MODEL="${TMP_BASE}_repeat.model"

#### - Pass2 - haploid
HAP_ONE_TMP="${TMP_BASE}_hap1_tmp.vcf.gz"
HAP_ONE_STD_TMP="${TMP_BASE}_hap1_nohp_tmp.vcf.gz"
HAP_ONE_PATCH="${TMP_BASE}_hap1_patch.vcf.gz"
HAP_ONE_OUT="${TMP_BASE}_hap1.vcf.gz"
HAP_TWO_TMP="${TMP_BASE}_hap2_tmp.vcf.gz"
HAP_TWO_STD_TMP="${TMP_BASE}_hap2_nohp_tmp.vcf.gz"
HAP_TWO_PATCH="${TMP_BASE}_hap2_patch.vcf.gz"
HAP_TWO_OUT="${TMP_BASE}_hap2.vcf.gz"

#### - Pass2 - diploid
DIPLOID_UNPHASED_HP="${TMP_BASE}_diploid_unphased_hp.vcf.gz"
DIPLOID_UNPHASED_PATCH="${TMP_BASE}_diploid_unphased_patch.vcf.gz"
DIPLOID_UNPHASED="${TMP_BASE}_diploid_unphased.vcf.gz"

#### - Merge generated VCF files
OUTPUT_GVCF="${_arg_output_vcf%%.vcf.gz}.g.vcf.gz"


#### Supporting scripts ####
DIR=$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd)
VCF_MOD_PY="${DIR}"/vcf_mod.py
GVCF_COMBINE="${DIR}"/gvcf_combine.py


set -e
#  Un-comment this line for more log information in the Bash shell
#set -exvuo pipefail


#### Pipeline fuctions ####
dnascope_hp()
{
    bam="$1" model="$2" repeat_model="$3" vcf="$4"
    bed="$5" read_filter="${6:-}" ds_model="${7:-}"
    ds_vcf="${8:-}"

    sentieon driver -t "$_arg_threads" -r "$_arg_reference_fasta" \
        --interval "$bed" -i "$bam" ${read_filter:+--read_filter $read_filter} \
        ${ds_vcf:+--algo DNAscope ${_arg_dbsnp:+--dbsnp "$_arg_dbsnp"} --model "$ds_model" "$ds_vcf"} \
        --algo DNAscopeHP ${_arg_dbsnp:+--dbsnp "$_arg_dbsnp"} \
        --model "$model" --pcr_indel_model "$repeat_model" \
        --min_repeat_count 6 "$vcf"
}

model_apply()
{
    input_vcf="$1"
    output_vcf="$2"
    model="$3"

    sentieon driver -t "$_arg_threads" -r "$_arg_reference_fasta" \
        --algo DNAModelApply --model "$model" -v "$input_vcf" "$output_vcf"
}

#### Pipeline ####
# Pass 1
## Call variants
sentieon driver -t "$_arg_threads" -r "$_arg_reference_fasta" \
    ${_arg_regions:+--interval "$_arg_regions"} -i "$_arg_sample_input" \
    ${_arg_gvcf:+--algo DNAscope --model "$_arg_model_bundle"/gvcf_model --emit_mode gvcf "$DIPLOID_GVCF"} \
    --algo DNAscope ${_arg_dbsnp:+--dbsnp "$_arg_dbsnp"} \
    --model "$_arg_model_bundle"/diploid_model "$DIPLOID_TMP_OUT"

model_apply "$DIPLOID_TMP_OUT" "$DIPLOID_OUT" "$_arg_model_bundle"/diploid_model

# Phasing
sentieon driver -t "$_arg_threads" -r "$_arg_reference_fasta" \
    -i "$_arg_sample_input" --algo VariantPhaser -v "$DIPLOID_OUT" \
    --max_depth 1000 --out_bed "$PHASED_BED" --out_ext "$PHASED_EXT" \
    "$PHASED_VCF"
if [ -n "$_arg_regions" ]; then
    bedtools subtract -a "$_arg_regions" -b "$PHASED_BED" > "$UNPHASED_BED" &
else
    (cat "$_arg_reference_fasta".fai | awk -v OFS='\t' '{print $1,0,$2}' > "$REF_BED"; \
        bedtools subtract -a "$REF_BED" -b "$PHASED_BED" > "$UNPHASED_BED") &
fi

## Create the repeat model
sentieon driver -t "$_arg_threads" -r "$_arg_reference_fasta" \
    -i "$_arg_sample_input" --interval "$PHASED_BED" \
    --read_filter PhasedReadFilter,phased_vcf="$PHASED_EXT",phase_select=tag \
    --algo RepeatModel --phased --min_map_qual 1 --min_group_count 10000 \
    --read_flag_mask drop=supplementary --repeat_extension 5 \
    --max_repeat_unit_size 2 --min_repeat_count 6 "$REPEAT_MODEL"

wait
bcftools view -T "$UNPHASED_BED" "$PHASED_VCF" | \
    sentieon util vcfconvert - "$PHASED_UNPHASED_VCF" &

# Pass 2 - call variants on the phased haploid chromosomes
## Call variants
dnascope_hp "$_arg_sample_input" "$_arg_model_bundle"/haploid_hp_model \
    "$REPEAT_MODEL" "$HAP_ONE_TMP" "$PHASED_BED" \
    "PhasedReadFilter,phased_vcf=$PHASED_EXT,phase_select=1" \
    "$_arg_model_bundle"/haploid_model "$HAP_ONE_STD_TMP"
dnascope_hp "$_arg_sample_input" "$_arg_model_bundle"/haploid_hp_model \
    "$REPEAT_MODEL" "$HAP_TWO_TMP" "$PHASED_BED" \
    "PhasedReadFilter,phased_vcf=$PHASED_EXT,phase_select=2" \
    "$_arg_model_bundle"/haploid_model "$HAP_TWO_STD_TMP"

## Merge DNAscope and DNAscopeHP VCFs
sentieon pyexec "$VCF_MOD_PY" -t "$_arg_threads" haploid_patch \
    --patch1 "$HAP_ONE_PATCH" --patch2 "$HAP_TWO_PATCH" \
    --hap1 "$HAP_ONE_STD_TMP" --hap2 "$HAP_TWO_STD_TMP" \
    --hap1_hp "$HAP_ONE_TMP" --hap2_hp "$HAP_TWO_TMP"

## Apply the trained model to the patched VCFs
model_apply "$HAP_ONE_PATCH" "$HAP_ONE_OUT" "$_arg_model_bundle"/haploid_model
model_apply "$HAP_TWO_PATCH" "$HAP_TWO_OUT" "$_arg_model_bundle"/haploid_model

# Pass 2 - call variant on the unphased regions
dnascope_hp "$_arg_sample_input" "$_arg_model_bundle"/diploid_hp_model \
    "$REPEAT_MODEL" "$DIPLOID_UNPHASED_HP" "$UNPHASED_BED"

## Merge the DNAscope and DNAscopeHP VCFs
wait
sentieon pyexec "$VCF_MOD_PY" -t "$_arg_threads" patch \
    --vcf "$PHASED_UNPHASED_VCF" --vcf_hp "$DIPLOID_UNPHASED_HP" \
    "$DIPLOID_UNPHASED_PATCH"
model_apply "$DIPLOID_UNPHASED_PATCH" "$DIPLOID_UNPHASED" \
    "$_arg_model_bundle"/diploid_model_unphased

## Merge the calls to create the output
sentieon pyexec "$VCF_MOD_PY" -t "$_arg_threads" merge \
    --hap1 "$HAP_ONE_OUT" --hap2 "$HAP_TWO_OUT" --unphased "$DIPLOID_UNPHASED" \
    --phased "$PHASED_VCF" --bed "$PHASED_BED" "${_arg_output_vcf}"

if [ "${_arg_gvcf}" = "on" ]; then
    sentieon pyexec "$GVCF_COMBINE" -t "$_arg_threads" "$DIPLOID_GVCF" "${_arg_output_vcf}" - | \
        sentieon util vcfconvert - "$OUTPUT_GVCF"
fi

rm -r "$tmp_dir"
