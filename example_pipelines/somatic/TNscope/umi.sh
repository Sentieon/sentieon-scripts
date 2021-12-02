#!/bin/sh

# Copyright (c) 2016-2020 Sentieon Inc. All rights reserved

# *******************************************
# Script to perform TNscope variant calling on
# UMI-tagged samples
# *******************************************

set -eu

# Update with the fullpath location of your sample fastq
TUMOR_SM="tumor_sample" #sample name
TUMOR_RGID="rg_$TUMOR_SM" #read group ID
PL="ILLUMINA" #or other sequencing platform
FASTQ_FOLDER="/home/pipeline/samples"
TUMOR_FASTQ_1="$FASTQ_FOLDER/tumor_1.fastq.gz"
TUMOR_FASTQ_2="$FASTQ_FOLDER/tumor_2.fastq.gz" #If using Illumina paired data

# Update with the location of the reference data files
FASTA_DIR="/home/regression/references/b37/"
FASTA="$FASTA_DIR/human_g1k_v37_decoy.fasta"
KNOWN_DBSNP="$FASTA_DIR/dbsnp_138.b37.vcf.gz"
INTERVAL_FILE="$FASTA_DIR/TruSeq_exome_targeted_regions.b37.bed"

# Update with the location of the Sentieon software package and license file
SENTIEON_INSTALL_DIR=/home/release/sentieon-genomics-|release_version|
export SENTIEON_LICENSE=/home/Licenses/Sentieon.lic #or using licsrvr: c1n11.sentieon.com:5443

#UMI information
READ_STRUCTURE="12M11S+T,+T" #an example duplex: "3M2S+T,3M2S+T" where duplex UMI extraction requires an identical read structure for both strands
DUPLEX_UMI="false" #set to "true" if duplex

# Other settings
NT=$(nproc) #number of threads to use in computation, set to number of cores in the server
START_DIR="$PWD/test/TNscope_umi" #Determine where the output files will be stored
BCFTOOLS_BINARY="/home/release/other_tools/bcftools-1.7/bcftools" #bcftools >=1.7 is recommended
MIN_DEPTH=100

# ******************************************
# 0. Setup
# ******************************************
WORKDIR="$START_DIR/$TUMOR_SM"
mkdir -p $WORKDIR
LOGFILE=$WORKDIR/run.log
exec >$LOGFILE 2>&1
cd $WORKDIR

# ******************************************
# 1. Pre-processing of FASTQ containing UMIs
# ******************************************
if [ "$DUPLEX_UMI" = "true" ] ; then
    READ_STRUCTURE="-d $READ_STRUCTURE"
fi
( $SENTIEON_INSTALL_DIR/bin/sentieon umi extract $READ_STRUCTURE $TUMOR_FASTQ_1 $TUMOR_FASTQ_2 || \
    { echo -n 'Extract error' >&2; exit -1; } ) | \
  ( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -p -C \
  -R "@RG\tID:$TUMOR_RGID\tSM:$TUMOR_SM\tPL:$PL" -t $NT \
  -K 10000000 $FASTA - || { echo -n 'BWA error'; exit 1; } ) | \
  $SENTIEON_INSTALL_DIR/bin/sentieon umi consensus -o umi_consensus.fastq.gz || \
  { echo "Alignment/Consensus failed"; exit 1; }

( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -p -C \
    -R "@RG\tID:$TUMOR_RGID\tSM:$TUMOR_SM\tPL:$PL" -t $NT -K 10000000 \
    $FASTA umi_consensus.fastq.gz || { echo -n 'BWA error'; exit 1; } ) | \
    $SENTIEON_INSTALL_DIR/bin/sentieon util sort --umi_post_process --sam2bam -i - \
    -o umi_consensus.bam || { echo "Consensus alignment failed"; exit 1; }

# ******************************************
# 2. Somatic and Structural variant calling
# ******************************************
# Consider adding `--disable_detector sv --trim_soft_clip` if not interested in SV calling
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i umi_consensus.bam \
    ${INTERVAL_FILE:+--interval_padding 10 --interval $INTERVAL_FILE} \
    --algo TNscope \
    --tumor_sample $TUMOR_SM \
    --dbsnp $KNOWN_DBSNP \
    --pcr_indel_model NONE \
    --disable_detector sv \
    --min_tumor_allele_frac 0.0003 \
    --filter_t_alt_frac 0.0003 \
    --max_error_per_read 3 \
    --min_init_tumor_lod 3.0 \
    --min_tumor_lod 3.0 \
    --min_base_qual 40 \
    --resample_depth 100000 \
    --assemble_mode 4 \
    output_tnscope.pre_filter.vcf.gz || \
    { echo "TNscope failed"; exit 1; }

#### The output of TNscope requires filtering to remove false positives.
#### Filter design depends on the specific sample and user needs to modify the following accordingly.
( $BCFTOOLS_BINARY annotate ${INTERVAL_FILE:+-R $INTERVAL_FILE} -x "FILTER/triallelic_site" output_tnscope.pre_filter.vcf.gz || \
        { echo "VCF filtering failed"; exit 1; } ) | \
    ( $BCFTOOLS_BINARY filter -m + -s "low_qual" -e "QUAL < 10" || \
        { echo "VCF filtering failed"; exit 1; } ) | \
    ( $BCFTOOLS_BINARY filter -m + -s "short_tandem_repeat" -e "RPA[0]>=10" || \
        { echo "VCF filtering failed"; exit 1; } ) | \
    ( $BCFTOOLS_BINARY filter -m + -s "read_pos_bias" -e "FMT/ReadPosRankSumPS[0] < -5" || \
        { echo "VCF filtering failed"; exit 1; } ) | \
    ( $BCFTOOLS_BINARY filter -m + -s "base_qual_bias" -e "FMT/BaseQRankSumPS[0] < -5" || \
        { echo "VCF filtering failed"; exit 1; } ) | \
    ( $BCFTOOLS_BINARY filter -m + -s "low_depth" -e "SUM(FMT/AD[0:]) < $MIN_DEPTH" || \
        { echo "VCF filtering failed"; exit 1; } ) | \
    $SENTIEON_INSTALL_DIR/bin/sentieon util vcfconvert - output_tnscope.filtered.vcf.gz || \
        { echo "VCF filtering failed"; exit 1; }
