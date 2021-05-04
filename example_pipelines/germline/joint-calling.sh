#!/bin/sh

# Copyright (c) 2016-2020 Sentieon Inc. All rights reserved

# *******************************************
# Script to perform joint calling on 3 samples
# with fastq files named sample<i>_1.fastq.gz
# and sample<i>_2.fastq.gz
# *******************************************

# Update with the fullpath location of your sample fastq
SM_LIST="sample1 sample2 sample3" # list of sample names
PL="ILLUMINA" #or other sequencing platform
FASTQ_FOLDER="/home/pipeline/samples"
FASTQ_1_SUFFIX="1.fastq.gz"
FASTQ_2_SUFFIX="2.fastq.gz" #If using Illumina paired data

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
NT=$(nproc) #number of threads to use in computation, set to number of cores in the server
START_DIR="$PWD/test/DNAseq_joint" #Determine where the output files will be stored





# You do not need to modify any of the lines below unless you want to tweak the pipeline

# ************************************************************************************************************************************************************************

# ******************************************
# 0. Setup
# ******************************************
WORKDIR="$START_DIR/joint_call"
mkdir -p $WORKDIR
LOGFILE=$WORKDIR/run.log
exec >$LOGFILE 2>&1

# ******************************************
# 0. Process all samples independently
# ******************************************
for SM in $SM_LIST; do
  RGID="rg_$SM"
  mkdir $WORKDIR/$SM
  cd $WORKDIR/$SM

  # ******************************************
  # 1. Mapping reads with BWA-MEM, sorting
  # ******************************************
  ( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -M -R "@RG\tID:$RGID\tSM:$SM\tPL:$PL" \
      -t $NT -K 10000000 $FASTA $FASTQ_FOLDER/${SM}_$FASTQ_1_SUFFIX \
      $FASTQ_FOLDER/${SM}_$FASTQ_2_SUFFIX || echo -n 'error' ) | \
      $SENTIEON_INSTALL_DIR/bin/sentieon util sort -r $FASTA -o sorted.bam -t $NT --sam2bam -i -
  if [ "$?" -ne "0" ]; then   echo "Alignment failed";   exit 1; fi
  
  # ******************************************
  # 2. Metrics
  # ******************************************
  $SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i sorted.bam \
      --algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution qd_metrics.txt \
      --algo GCBias --summary gc_summary.txt gc_metrics.txt --algo AlignmentStat \
      --adapter_seq '' aln_metrics.txt --algo InsertSizeMetricAlgo is_metrics.txt
  if [ "$?" -ne "0" ]; then   echo "Metrics failed";   exit 1; fi

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
  if [ "$?" -ne "0" ]; then   echo "LocusCollector failed";   exit 1; fi

  $SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NT -i sorted.bam --algo Dedup \
      --score_info score.txt --metrics dedup_metrics.txt deduped.bam
  if [ "$?" -ne "0" ]; then   echo "Dedup failed";   exit 1; fi
  
  # ******************************************
  # 2a. Coverage metrics
  # ******************************************
  $SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i deduped.bam \
      --algo CoverageMetrics coverage_metrics
  if [ "$?" -ne "0" ]; then   echo "CoverageMetrics failed";   exit 1; fi

  # ******************************************
  # 5. Base recalibration
  # ******************************************
  $SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i deduped.bam --algo QualCal \
      -k $KNOWN_DBSNP -k $KNOWN_MILLS -k $KNOWN_INDELS recal_data.table
  if [ "$?" -ne "0" ]; then   echo "QualCal1 failed";   exit 1; fi

  $SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i deduped.bam -q recal_data.table \
      --algo QualCal -k $KNOWN_DBSNP -k $KNOWN_MILLS -k $KNOWN_INDELS recal_data.table.post
  if [ "$?" -ne "0" ]; then   echo "QualCal2 failed";   exit 1; fi

  $SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NT --algo QualCal --plot \
      --before recal_data.table --after recal_data.table.post recal.csv
  if [ "$?" -ne "0" ]; then   echo "QualCal3 failed";   exit 1; fi

  $SENTIEON_INSTALL_DIR/bin/sentieon plot QualCal -o recal_plots.pdf recal.csv
  if [ "$?" -ne "0" ]; then   echo "QualCal4 failed";   exit 1; fi
  
  # ******************************************
  # 6b. HC Variant caller for GVCF
  # ******************************************
  $SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i deduped.bam -q recal_data.table \
      --algo Haplotyper -d $KNOWN_DBSNP --emit_conf=30 --call_conf=30 --emit_mode gvcf \
      output-hc.g.vcf.gz
  if [ "$?" -ne "0" ]; then   echo "Haplotyper failed";   exit 1; fi

  GVCF_INPUTS="$GVCF_INPUTS -v $WORKDIR/$SM/output-hc.g.vcf.gz"
done

# ******************************************
# Perform the joint calling 
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA --algo GVCFtyper $GVCF_INPUTS \
    -d $KNOWN_DBSNP $WORKDIR/output-jc.vcf.gz
if [ "$?" -ne "0" ]; then   echo "GVCFtyper failed";   exit 1; fi
