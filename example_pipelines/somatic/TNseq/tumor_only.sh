#!/bin/sh

# Copyright (c) 2016-2020 Sentieon Inc. All rights reserved

# *******************************************
# Script to perform TN seq variant calling
# using a Tumor sample with fastq files
# named tumor_1.fastq.gz, tumor_2.fastq.gz
# and data from a Panel of Normals and Cosmic DB
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
KNOWN_INDELS="$FASTA_DIR/1000G_phase1.indels.b37.vcf.gz"
KNOWN_MILLS="$FASTA_DIR/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
# Update with the location of the panel of normal and CosmicDB vcf files
# We recommend that you create the panel of normal file with the corresponding algorithm that you plan to use for the somatic mutation calling. 
PANEL_OF_NORMAL_TNSNV=
PANEL_OF_NORMAL_TNHAPLOTYPER=
PANEL_OF_NORMAL_TNHAPLOTYPER2=
COSMIC_DB="/home/regression/references/b37/b37_cosmic_v54_120711.vcf.gz"
CONTAMINATION_VCF="$FASTA_DIR/germline_vcf-af-only-gnomad.raw.sites.vcf"  # A VCF of germline sites to use for contamination detection
GERMLINE_VCF=  # A VCF of known germline sites

# Update with the location of the Sentieon software package and license file
SENTIEON_INSTALL_DIR=/home/release/sentieon-genomics-|release_version|
export SENTIEON_LICENSE=/home/Licenses/Sentieon.lic #or using licsrvr: c1n11.sentieon.com:5443

# Other settings
NT=$(nproc) #number of threads to use in computation, set to number of cores in the server
START_DIR="$PWD/test/TNseq_tumoronly" #Determine where the output files will be stored





# You do not need to modify any of the lines below unless you want to tweak the pipeline

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
    { echo "Alignment failed"; exit 1; }

# ******************************************
# 2a. Metrics for tumor sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i tumor_sorted.bam \
    --algo MeanQualityByCycle tumor_mq_metrics.txt \
    --algo QualDistribution tumor_qd_metrics.txt --algo GCBias \
    --summary tumor_gc_summary.txt tumor_gc_metrics.txt --algo AlignmentStat \
    --adapter_seq '' tumor_aln_metrics.txt \
    --algo InsertSizeMetricAlgo tumor_is_metrics.txt \
    --algo CoverageMetrics --omit_base_output tumor_coverage_metrics || \
    { echo "Metrics failed"; exit 1; }

$SENTIEON_INSTALL_DIR/bin/sentieon plot GCBias -o tumor_gc-report.pdf tumor_gc_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualDistribution \
    -o tumor_qd-report.pdf tumor_qd_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot MeanQualityByCycle \
    -o tumor_mq-report.pdf tumor_mq_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot InsertSizeMetricAlgo \
    -o tumor_is-report.pdf tumor_is_metrics.txt

# ******************************************
# 3a. Remove Duplicate Reads for tumor
# sample. It is possible
# to mark instead of remove duplicates
# by ommiting the --rmdup option in Dedup
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NT -i tumor_sorted.bam --algo LocusCollector \
    --fun score_info tumor_score.txt || { echo "LocusCollector failed"; exit 1; }

$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NT -i tumor_sorted.bam --algo Dedup \
    --rmdup --score_info tumor_score.txt --metrics tumor_dedup_metrics.txt tumor_deduped.bam || \
    { echo "Dedup failed"; exit 1; }

# ******************************************
# 4a. Somatic Variant Calling - TNhaplotyper2
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i tumor_deduped.bam \
    --algo TNhaplotyper2 --tumor_sample $TUMOR_SM \
    ${PANEL_OF_NORMAL_TNHAPLOTYPER2:+--pon $PANEL_OF_NORMAL_TNHAPLOTYPER2} \
    ${GERMLINE_VCF:+--germline_vcf $GERMLINE_VCF} output-tnhap2-tmp.vcf.gz \
    --algo OrientationBias --tumor_sample $TUMOR_SM output-orientation \
    --algo ContaminationModel --tumor_sample $TUMOR_SM  --vcf $CONTAMINATION_VCF \
    --tumor_segments output-contamination-segments output-contamination || \
    { echo "TNhaplotyper2 failed"; exit 1; }

$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA --algo TNfilter \
    -v output-tnhap2-tmp.vcf.gz --tumor_sample $TUMOR_SM \
    --contamination output-contamination --tumor_segments output-contamination-segments \
    --orientation_priors output-orientation output-tnhap2.vcf.gz || \
    { echo "TNfilter failed"; exit 1; }


# Uncomment the following commands to run somatic variant calling with TNhaplotyper
# ******************************************
# 4b. Somatic Variant Calling - TNhaplotyper
# ******************************************
#$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i tumor_deduped.bam \
#    --algo TNhaplotyper --tumor_sample $TUMOR_SM \
#    ${PANEL_OF_NORMAL_TNHAPLOTYPER:+--pon $PANEL_OF_NORMAL_TNHAPLOTYPER} \
#    --cosmic $COSMIC_DB --dbsnp $KNOWN_DBSNP \
#    output-tnhaplotyper.vcf.gz || { echo "TNhaplotyper failed"; exit 1; }


# Uncomment the following commands to run indel realignment, and somatic variant
# calling with TNsnv
# ******************************************
# 5. Indel realigner for tumor sample
# ******************************************
#$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i tumor_deduped.bam \
#    --algo Realigner -k $KNOWN_MILLS -k $KNOWN_INDELS tumor_realigned.bam || \
#    { echo "Realigner failed"; exit 1; }

# ******************************************
# 6. Somatic Variant Calling - TNsnv
# ******************************************
#$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i tumor_realigned.bam \
#    --algo TNsnv --tumor_sample $TUMOR_SM \
#    ${PANEL_OF_NORMAL_TNSNV:+--pon $PANEL_OF_NORMAL_TNSNV} \
#    --cosmic $COSMIC_DB --dbsnp $KNOWN_DBSNP \
#    --call_stats_out output-call.stats output-tnsnv.vcf.gz || \
#    { echo "TNsnv failed"; exit 1; }

