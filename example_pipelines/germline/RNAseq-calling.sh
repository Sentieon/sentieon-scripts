#!/bin/sh

# Copyright (c) 2016-2020 Sentieon Inc. All rights reserved

# *******************************************
# Script to perform RNA variant calling
# using a single sample with fastq files
# named 1.fastq.gz and 2.fastq.gz
# *******************************************

# Update with the fullpath location of your sample fastq
SM="sample" #sample name
RGID="rg_$SM" #read group ID
PL="ILLUMINA" #or other sequencing platform
FASTQ_FOLDER="/home/pipeline/samples"
FASTQ_1="$FASTQ_FOLDER/rna1.fastq.gz"
FASTQ_2="$FASTQ_FOLDER/rna2.fastq.gz" #If using Illumina paired data

# Update with the location of the reference data files
FASTA_DIR="/home/regression/references/b37/"
FASTA="$FASTA_DIR/human_g1k_v37_decoy.fasta"
KNOWN_DBSNP="$FASTA_DIR/dbsnp_138.b37.vcf.gz"
KNOWN_INDELS="$FASTA_DIR/1000G_phase1.indels.b37.vcf.gz"
KNOWN_MILLS="$FASTA_DIR/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
#comment if STAR should generate a new genomeDir
STAR_FASTA="$FASTA_DIR/human_g1k_v37_decoy.fasta.genomeDir"
#uncomment if you would like to use a bed file
#INTERVAL_FILE="$FASTA_DIR/TruSeq_exome_targeted_regions.b37.bed"

# Update with the location of the Sentieon software package and license file
SENTIEON_INSTALL_DIR=/home/release/sentieon-genomics-|release_version|
export SENTIEON_LICENSE=/home/Licenses/Sentieon.lic #or using licsrvr: c1n11.sentieon.com:5443

# Other settings
NT=$(nproc) #number of threads to use in computation, set to number of cores in the server
START_DIR="$PWD/test/RNAseq" #Determine where the output files will be stored
STAR_BINARY="/home/release/other_tools/STAR/bin/Linux_x86_64_static/STAR"





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
DRIVER_INTERVAL_OPTION="${INTERVAL_FILE:+--interval $INTERVAL_FILE}"

# ******************************************
# 1. Mapping reads with STAR
# ******************************************
if [ -z "$STAR_FASTA" ]; then
  STAR_FASTA="genomeDir"
  # The genomeDir generation could be reused
  mkdir $star_fasta
  $STAR_BINARY --runMode genomeGenerate --genomeDir $STAR_FASTA --genomeFastaFiles $FASTA \
      --runThreadN $NT
fi
#perform the actual alignment and sorting
$STAR_BINARY --twopassMode Basic --genomeDir $STAR_FASTA --runThreadN $NT --outSAMtype BAM \
    SortedByCoordinate --twopass1readsN -1 --sjdbOverhang 75 --readFilesIn $FASTQ_1 $FASTQ_2 \
    --readFilesCommand zcat --outSAMattrRGline ID:$RGID SM:$SM PL:$PL
mv Aligned.sortedByCoord.out.bam sorted.bam
$SENTIEON_INSTALL_DIR/bin/sentieon util index sorted.bam

# ******************************************
# 2. Metrics
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver $DRIVER_INTERVAL_OPTION -r $FASTA -t $NT \
    -i sorted.bam --algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution \
    qd_metrics.txt --algo GCBias --summary gc_summary.txt gc_metrics.txt \
    --algo AlignmentStat --adapter_seq '' aln_metrics.txt --algo InsertSizeMetricAlgo \
    is_metrics.txt
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
    --fun score_info score.txt
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NT -i sorted.bam --algo Dedup \
    --score_info score.txt --metrics dedup_metrics.txt deduped.bam

# ******************************************
# 2a. Coverage metrics
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i deduped.bam \
    --algo CoverageMetrics coverage_metrics

# ******************************************
# 4. Split reads at Junction
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i deduped.bam \
    --algo RNASplitReadsAtJunction --reassign_mapq 255:60 splitted.bam

# ******************************************
# 6. Base recalibration
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver $DRIVER_INTERVAL_OPTION -r $FASTA -t $NT \
    -i splitted.bam --algo QualCal -k $KNOWN_DBSNP -k $KNOWN_MILLS -k $KNOWN_INDELS \
    recal_data.table
$SENTIEON_INSTALL_DIR/bin/sentieon driver $DRIVER_INTERVAL_OPTION -r $FASTA -t $NT \
    -i splitted.bam -q recal_data.table --algo QualCal -k $KNOWN_DBSNP -k $KNOWN_MILLS \
    -k $KNOWN_INDELS recal_data.table.post
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NT --algo QualCal --plot \
    --before recal_data.table --after recal_data.table.post recal.csv
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualCal -o recal_plots.pdf recal.csv

# ******************************************
# 7. HC Variant caller for RNA
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver $DRIVER_INTERVAL_OPTION -r $FASTA -t $NT \
    -i splitted.bam -q recal_data.table --algo Haplotyper -d $KNOWN_DBSNP \
    --trim_soft_clip --emit_conf=20 --call_conf=20 output-hc-rna.vcf.gz
