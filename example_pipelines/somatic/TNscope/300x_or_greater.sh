#!/bin/sh

# Copyright (c) 2016-2020 Sentieon Inc. All rights reserved

# *******************************************
# Script to perform TNscope variant calling on
# samples with coverage over 300x using a matched
# paired Tumor+normal sample with fastq files
# named normal_1.fastq.gz, normal_2.fastq.gz,
# tumor_1.fastq.gz, tumor_2.fastq.gz
# *******************************************

set -eu

# Update with the fullpath location of your sample fastq
TUMOR_SM="tumor_sample" #sample name
TUMOR_RGID="rg_$TUMOR_SM" #read group ID
NORMAL_SM="normal_sample" #sample name
NORMAL_RGID="rg_$NORMAL_SM" #read group ID
PL="ILLUMINA" #or other sequencing platform
FASTQ_FOLDER="/home/pipeline/samples"
TUMOR_FASTQ_1="$FASTQ_FOLDER/tumor_1.fastq.gz"
TUMOR_FASTQ_2="$FASTQ_FOLDER/tumor_2.fastq.gz" #If using Illumina paired data
NORMAL_FASTQ_1="$FASTQ_FOLDER/normal_1.fastq.gz"
NORMAL_FASTQ_2="$FASTQ_FOLDER/normal_2.fastq.gz" #If using Illumina paired data

# Update with the location of the reference data files
FASTA_DIR="/home/regression/references/b37/"
FASTA="$FASTA_DIR/human_g1k_v37_decoy.fasta"
KNOWN_DBSNP="$FASTA_DIR/dbsnp_138.b37.vcf.gz"
KNOWN_INDELS="$FASTA_DIR/1000G_phase1.indels.b37.vcf.gz"
KNOWN_MILLS="$FASTA_DIR/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
INTERVAL_FILE="$FASTA_DIR/TruSeq_exome_targeted_regions.b37.bed"

# Update with the location of the Sentieon software package and license file
SENTIEON_INSTALL_DIR=/home/release/sentieon-genomics-|release_version|
export SENTIEON_LICENSE=/home/Licenses/Sentieon.lic #or using licsrvr: c1n11.sentieon.com:5443

# Other settings
NT=$(nproc) #number of threads to use in computation, set to number of cores in the server
START_DIR="$PWD/test/TNscope_highcov" #Determine where the output files will be stored
BCFTOOLS_BINARY="/home/release/other_tools/bcftools-1.7/bcftools" #bcftools >=1.7 is recommended
MIN_QUAL=10







# ************************************************************************************************************************************************************************

# ******************************************
# 0. Setup
# ******************************************
WORKDIR="$START_DIR"
mkdir -p $WORKDIR
LOGFILE=$WORKDIR/run.log
exec >$LOGFILE 2>&1
cd $WORKDIR

# ******************************************
# 1a. Mapping reads with BWA-MEM, sorting for tumor sample
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000 
( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -R "@RG\tID:$TUMOR_RGID\tSM:$TUMOR_SM\tPL:$PL" \
    -t $NT -K 10000000 $FASTA $TUMOR_FASTQ_1 $TUMOR_FASTQ_2 || \
    { echo -n 'BWA error'; exit 1; } ) | \
    $SENTIEON_INSTALL_DIR/bin/sentieon util sort -o tumor_sorted.bam -t $NT --sam2bam -i - || \
    { echo "Alignment1 failed"; exit 1; }

# ******************************************
# 1b. Mapping reads with BWA-MEM, sorting for normal sample
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000 
( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -R "@RG\tID:$NORMAL_RGID\tSM:$NORMAL_SM\tPL:$PL" \
    -t $NT -K 10000000 $FASTA $NORMAL_FASTQ_1 $NORMAL_FASTQ_2 || \
    { echo -n 'BWA error'; exit 1; } ) | \
    $SENTIEON_INSTALL_DIR/bin/sentieon util sort -o normal_sorted.bam -t $NT --sam2bam -i - || \
    { echo "Alignment2 failed"; exit 1; }

# ******************************************
# 2a. Metrics for tumor sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i tumor_sorted.bam \
    --algo MeanQualityByCycle tumor_mq_metrics.txt \
    --algo QualDistribution tumor_qd_metrics.txt --algo GCBias \
    --summary tumor_gc_summary.txt tumor_gc_metrics.txt --algo AlignmentStat \
    --adapter_seq '' tumor_aln_metrics.txt --algo InsertSizeMetricAlgo tumor_is_metrics.txt || \
    { echo "Metrics1 failed"; exit 1; }

$SENTIEON_INSTALL_DIR/bin/sentieon plot GCBias -o tumor_gc-report.pdf tumor_gc_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualDistribution \
    -o tumor_qd-report.pdf tumor_qd_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot MeanQualityByCycle \
    -o tumor_mq-report.pdf tumor_mq_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot InsertSizeMetricAlgo \
    -o tumor_is-report.pdf tumor_is_metrics.txt

# ******************************************
# 2b. Metrics for normal sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i normal_sorted.bam \
    --algo MeanQualityByCycle normal_mq_metrics.txt \
    --algo QualDistribution normal_qd_metrics.txt --algo GCBias \
    --summary normal_gc_summary.txt normal_gc_metrics.txt --algo AlignmentStat \
    --adapter_seq '' normal_aln_metrics.txt --algo InsertSizeMetricAlgo normal_is_metrics.txt || \
    { echo "Metrics2 failed"; exit 1; }

$SENTIEON_INSTALL_DIR/bin/sentieon plot GCBias -o normal_gc-report.pdf normal_gc_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualDistribution \
    -o normal_qd-report.pdf normal_qd_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot MeanQualityByCycle \
    -o normal_mq-report.pdf normal_mq_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot InsertSizeMetricAlgo \
    -o normal_is-report.pdf normal_is_metrics.txt

# ******************************************
# 3a. Remove Duplicate Reads for tumor
# sample. It is possible
# to remove instead of mark duplicates
# by adding the --rmdup option in Dedup
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NT -i tumor_sorted.bam --algo LocusCollector \
    --fun score_info tumor_score.txt || { echo "LocusCollector1 failed"; exit 1; }

$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NT -i tumor_sorted.bam --algo Dedup \
    --score_info tumor_score.txt --metrics tumor_dedup_metrics.txt tumor_deduped.bam || \
    { echo "Dedup1 failed"; exit 1; }

# ******************************************
# 3b. Remove Duplicate Reads for normal
# sample. It is possible
# to remove instead of mark duplicates
# by adding the --rmdup option in Dedup
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NT -i normal_sorted.bam --algo LocusCollector \
    --fun score_info normal_score.txt || { echo "LocusCollector2 failed"; exit 1; }

$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NT -i normal_sorted.bam --algo Dedup \
    --score_info normal_score.txt --metrics normal_dedup_metrics.txt normal_deduped.bam || \
    { echo "Dedup2 failed"; exit 1; }

# ******************************************
# 2a. Coverage metrics for normal sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i normal_deduped.bam \
    --algo CoverageMetrics normal_coverage_metrics || { echo "CoverageMetrics1 failed"; exit 1; }

# ******************************************
# 2b. Coverage metrics for tumor sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i tumor_deduped.bam \
    --algo CoverageMetrics tumor_coverage_metrics || { echo "CoverageMetrics2 failed"; exit 1; }

# ******************************************
# 4a. Base recalibration for tumor sample (Skip if Panel or targeted sequencing)
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i tumor_deduped.bam \
    --algo QualCal -k $KNOWN_DBSNP -k $KNOWN_MILLS -k $KNOWN_INDELS tumor_recal_data.table
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i tumor_deduped.bam \
    -q tumor_recal_data.table --algo QualCal -k $KNOWN_DBSNP -k $KNOWN_MILLS \
    -k $KNOWN_INDELS tumor_recal_data.table.post
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NT --algo QualCal --plot \
    --before tumor_recal_data.table --after tumor_recal_data.table.post tumor_recal.csv
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualCal -o tumor_recal_plots.pdf tumor_recal.csv

# ******************************************
# 4b. Base recalibration for normal sample (Skip if Panel or targeted sequencing)
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i normal_deduped.bam \
    --algo QualCal -k $KNOWN_DBSNP -k $KNOWN_MILLS -k $KNOWN_INDELS normal_recal_data.table
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i normal_deduped.bam \
    -q normal_recal_data.table --algo QualCal -k $KNOWN_DBSNP -k $KNOWN_MILLS \
    -k $KNOWN_INDELS normal_recal_data.table.post
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NT --algo QualCal --plot \
    --before normal_recal_data.table --after normal_recal_data.table.post normal_recal.csv
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualCal -o normal_recal_plots.pdf normal_recal.csv

# ******************************************
# 5. Somatic and Structural variant calling
# ******************************************
# Consider adding `--disable_detector sv --trim_soft_clip` if not interested in SV calling
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT \
    ${INTERVAL_FILE:+--interval_padding 250 --interval $INTERVAL_FILE} \
    -i tumor_deduped.bam -i normal_deduped.bam \
    -q tumor_recal_data.table -q normal_recal_data.table \
    --algo TNscope \
    --prune_factor 0 \
    --tumor_sample $TUMOR_SM --normal_sample $NORMAL_SM --dbsnp $KNOWN_DBSNP \
    --sv_mask_ext 10 --max_fisher_pv_active 0.05 --min_tumor_allele_frac 0.005 \
    --filter_t_alt_frac 0.005 --resample_depth 100000 --clip_by_minbq 1 \
    --assemble_mode 4 output_tnscope.pre_filter.vcf.gz || \
    { echo "TNscope failed"; exit 1; }

# ******************************************
# 6. Variant filtration
# ******************************************
( $BCFTOOLS_BINARY annotate ${INTERVAL_FILE:+-R $INTERVAL_FILE} -x "FILTER/triallelic_site" output_tnscope.pre_filter.vcf.gz || \
        { echo "VCF filtering failed"; exit 1; } ) | \
    ( $BCFTOOLS_BINARY filter -m + -s "insignificant" -e "(PV>0.25 && PV2>0.25) || (INFO/STR == 1 && PV>0.05)" || \
        { echo "VCF filtering failed"; exit 1; } ) | \
    ( $BCFTOOLS_BINARY filter -m + -s "low_qual" -e "QUAL < $MIN_QUAL" || \
        { echo "VCF filtering failed"; exit 1; } ) | \
    ( $BCFTOOLS_BINARY filter -m + -s "short_tandem_repeat" -e "RPA[0]>=10" || \
        { echo "VCF filtering failed"; exit 1; } ) | \
    ( $BCFTOOLS_BINARY filter -m + -s "read_pos_bias" -e "FMT/ReadPosRankSumPS[0] < -5" || \
        { echo "VCF filtering failed"; exit 1; } ) | \
	$SENTIEON_INSTALL_DIR/bin/sentieon util vcfconvert - output_tnscope.filtered.vcf.gz || \
        { echo "VCF filtering failed"; exit 1; }
