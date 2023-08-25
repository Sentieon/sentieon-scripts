#!/bin/sh

# Copyright (c) 2016-2023 Sentieon Inc. All rights reserved

# *******************************************
# Script to perform DNAscope variant calling
# using a single sample with fastq files
# named 1.fastq.gz and 2.fastq.gz
# *******************************************

set -eu

# Update with the fullpath location of your sample fastq
SM="sample" #sample name
RGID="rg_$SM" #read group ID
PL="ILLUMINA" #or other sequencing platform
FASTQ_FOLDER="/home/pipeline/samples"
FASTQ_1="$FASTQ_FOLDER/1.fastq.gz"
FASTQ_2="$FASTQ_FOLDER/2.fastq.gz" #If using Illumina paired data

# Update with the location to the DNAscope model file
# Model files can be found at, https://github.com/Sentieon/sentieon-models
DNASCOPE_MODEL=/home/pipeline/models/DNAscopeIlluminaWGS2.0.bundle

# Update with the location of the reference data files
FASTA_DIR="/home/regression/references/b37/"
FASTA="$FASTA_DIR/human_g1k_v37_decoy.fasta"
KNOWN_DBSNP="$FASTA_DIR/dbsnp_138.b37.vcf.gz"
KNOWN_INDELS="$FASTA_DIR/1000G_phase1.indels.b37.vcf.gz"
KNOWN_MILLS="$FASTA_DIR/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"

# Update with the location of the Sentieon software package and license file
SENTIEON_INSTALL_DIR=/home/release/sentieon-genomics-|release_version|
export SENTIEON_LICENSE=/home/Licenses/Sentieon.lic #or using licsrvr: c1n11.sentieon.com:5443

# Other settings
PCRFREE=true # The data was sequenced with a PCR-free library prep
NT=$(nproc) #number of threads to use in computation, set to number of cores in the server
START_DIR="$PWD/test/DNAscope" #Determine where the output files will be stored





# You do not need to modify any of the lines below unless you want to tweak the pipeline

# ************************************************************************************************************************************************************************

# ******************************************
# 0. Setup
# ******************************************
WORKDIR="$START_DIR/${SM}" 
mkdir -p $WORKDIR
LOGFILE=$WORKDIR/run.log
exec >$LOGFILE 2>&1
cd $WORKDIR

# ******************************************
# 1. Mapping reads with BWA-MEM, sorting
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000 
( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -R "@RG\tID:$RGID\tSM:$SM\tPL:$PL" \
    -t $NT -K 10000000 -x $DNASCOPE_MODEL/bwa.model \
    $FASTA $FASTQ_1 $FASTQ_2 || { echo -n 'BWA error'; exit 1; } ) | \
    $SENTIEON_INSTALL_DIR/bin/sentieon util sort -r $FASTA -o sorted.bam -t $NT \
    --sam2bam -i - || { echo "Alignment failed"; exit 1; }

# ******************************************
# 2. Metrics
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i sorted.bam \
    --algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution qd_metrics.txt \
    --algo GCBias --summary gc_summary.txt gc_metrics.txt --algo AlignmentStat \
    --adapter_seq '' aln_metrics.txt --algo InsertSizeMetricAlgo is_metrics.txt || \
    { echo "Metrics failed"; exit 1; }

$SENTIEON_INSTALL_DIR/bin/sentieon plot GCBias -o gc-report.pdf gc_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualDistribution -o qd-report.pdf qd_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot MeanQualityByCycle -o mq-report.pdf mq_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot InsertSizeMetricAlgo -o is-report.pdf is_metrics.txt

# ******************************************
# 3. Remove Duplicate Reads. It is possible
# to remove instead of mark duplicates
# by adding the --rmdup option in Dedup
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NT -i sorted.bam --algo LocusCollector \
    --fun score_info score.txt || { echo "LocusCollector failed"; exit 1; }

$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NT -i sorted.bam --algo Dedup \
    --score_info score.txt --metrics dedup_metrics.txt deduped.bam || \
    { echo "Dedup failed"; exit 1; }

# ******************************************
# 4a. DNAscope variant calling
# ******************************************
indel_model_arg=""
if [ "$PCRFREE" = true ]; then
    indel_model_arg="--pcr_indel_model none"
fi
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i deduped.bam \
    --algo DNAscope $indel_model_arg --model $DNASCOPE_MODEL/dnascope.model \
    -d $KNOWN_DBSNP output-ds_tmp.vcf.gz || \
    { echo "DNAscope failed"; exit 1; }

# ******************************************
# 4b. Variant filtering and genotyping
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT \
    --algo DNAModelApply --model $DNASCOPE_MODEL/dnascope.model \
    -v output-ds_tmp.vcf.gz output-ds.vcf.gz || \
    { echo "DNAModelApply failed"; exit 1; }
