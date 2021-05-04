#!/bin/sh

# Copyright (c) 2016-2020 Sentieon Inc. All rights reserved

# *******************************************
# Script to perform TN seq variant calling
# using a matched paired Tumor+normal sample with fastq
# files named normal_1.fastq.gz, normal_2.fastq.gz
# tumor_1.fastq.gz, tumor_2.fastq.gz
# *******************************************


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
NORMAL_FASTQ_2="$FASTQ_FOLDER/normal_2.fastq.gz"

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
START_DIR="$PWD/test/TNseq" #Determine where the output files will be stored





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
( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -M -R "@RG\tID:$TUMOR_RGID\tSM:$TUMOR_SM\tPL:$PL" \
    -t $NT -K 10000000 $FASTA $TUMOR_FASTQ_1 $TUMOR_FASTQ_2 || echo -n 'error' ) | \
    $SENTIEON_INSTALL_DIR/bin/sentieon util sort -o tumor_sorted.bam -t $NT --sam2bam -i -
if [ "$?" -ne "0" ]; then   echo "Alignment1 failed";   exit 1; fi

# ******************************************
# 1b. Mapping reads with BWA-MEM, sorting for normal sample
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000 
( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -M -R "@RG\tID:$NORMAL_RGID\tSM:$NORMAL_SM\tPL:$PL" \
    -t $NT -K 10000000 $FASTA $NORMAL_FASTQ_1 $NORMAL_FASTQ_2 || echo -n 'error' ) | \
    $SENTIEON_INSTALL_DIR/bin/sentieon util sort -o normal_sorted.bam -t $NT --sam2bam -i -
if [ "$?" -ne "0" ]; then   echo "Alignment2 failed";   exit 1; fi

# ******************************************
# 2a. Metrics for tumor sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i tumor_sorted.bam \
    --algo MeanQualityByCycle tumor_mq_metrics.txt \
    --algo QualDistribution tumor_qd_metrics.txt --algo GCBias \
    --summary tumor_gc_summary.txt tumor_gc_metrics.txt --algo AlignmentStat \
    --adapter_seq '' tumor_aln_metrics.txt --algo InsertSizeMetricAlgo tumor_is_metrics.txt
if [ "$?" -ne "0" ]; then   echo "Metrics1 failed";   exit 1; fi

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
    --adapter_seq '' normal_aln_metrics.txt --algo InsertSizeMetricAlgo normal_is_metrics.txt
if [ "$?" -ne "0" ]; then   echo "Metrics2 failed";   exit 1; fi

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
    --fun score_info tumor_score.txt
if [ "$?" -ne "0" ]; then   echo "LocusCollector1 failed";   exit 1; fi

$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NT -i tumor_sorted.bam --algo Dedup \
    --score_info tumor_score.txt --metrics tumor_dedup_metrics.txt tumor_deduped.bam
if [ "$?" -ne "0" ]; then   echo "Dedup1 failed";   exit 1; fi

# ******************************************
# 3b. Remove Duplicate Reads for normal
# sample. It is possible
# to remove instead of mark duplicates
# by adding the --rmdup option in Dedup
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NT -i normal_sorted.bam --algo LocusCollector --fun score_info normal_score.txt
if [ "$?" -ne "0" ]; then   echo "LocusCollector2 failed";   exit 1; fi

$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NT -i normal_sorted.bam --algo Dedup --score_info normal_score.txt --metrics normal_dedup_metrics.txt normal_deduped.bam 
if [ "$?" -ne "0" ]; then   echo "Dedup2 failed";   exit 1; fi

# ******************************************
# 2a. Coverage metrics for tumor sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i tumor_deduped.bam \
    --algo CoverageMetrics tumor_coverage_metrics
if [ "$?" -ne "0" ]; then   echo "CoverageMetrics1 failed";   exit 1; fi

# ******************************************
# 2b. Coverage metrics for normal sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i normal_deduped.bam \
    --algo CoverageMetrics normal_coverage_metrics
if [ "$?" -ne "0" ]; then   echo "CoverageMetrics2 failed";   exit 1; fi

# ******************************************
# 4a. Indel realigner for tumor sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i tumor_deduped.bam \
    --algo Realigner -k $KNOWN_MILLS -k $KNOWN_INDELS tumor_realigned.bam
if [ "$?" -ne "0" ]; then   echo "Realigner1 failed";   exit 1; fi

# ******************************************
# 4b. Indel realigner for normal sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i normal_deduped.bam \
    --algo Realigner -k $KNOWN_MILLS -k $KNOWN_INDELS normal_realigned.bam
if [ "$?" -ne "0" ]; then   echo "Realigner2 failed";   exit 1; fi

# ******************************************
# 5a. Base recalibration for tumor sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i tumor_realigned.bam \
    --algo QualCal -k $KNOWN_DBSNP -k $KNOWN_MILLS -k $KNOWN_INDELS tumor_recal_data.table
if [ "$?" -ne "0" ]; then   echo "QualCal1 failed";   exit 1; fi

$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i tumor_realigned.bam \
    -q tumor_recal_data.table --algo QualCal -k $KNOWN_DBSNP -k $KNOWN_MILLS \
    -k $KNOWN_INDELS tumor_recal_data.table.post
if [ "$?" -ne "0" ]; then   echo "QualCal2 failed";   exit 1; fi

$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NT --algo QualCal --plot \
    --before tumor_recal_data.table --after tumor_recal_data.table.post tumor_recal.csv
if [ "$?" -ne "0" ]; then   echo "QualCal3 failed";   exit 1; fi

$SENTIEON_INSTALL_DIR/bin/sentieon plot QualCal -o tumor_recal_plots.pdf tumor_recal.csv
# ******************************************
# 5b. Base recalibration for normal sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i normal_realigned.bam \
    --algo QualCal -k $KNOWN_DBSNP -k $KNOWN_MILLS -k $KNOWN_INDELS normal_recal_data.table
if [ "$?" -ne "0" ]; then   echo "QualCal4 failed";   exit 1; fi

$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i normal_realigned.bam \
    -q normal_recal_data.table --algo QualCal -k $KNOWN_DBSNP -k $KNOWN_MILLS \
    -k $KNOWN_INDELS normal_recal_data.table.post
if [ "$?" -ne "0" ]; then   echo "QualCal5 failed";   exit 1; fi

$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NT --algo QualCal --plot \
    --before normal_recal_data.table --after normal_recal_data.table.post normal_recal.csv
if [ "$?" -ne "0" ]; then   echo "QualCal6 failed";   exit 1; fi

$SENTIEON_INSTALL_DIR/bin/sentieon plot QualCal -o normal_recal_plots.pdf normal_recal.csv

# ******************************************
# 6. Corealignment of tumor and normal
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i tumor_realigned.bam \
    -i normal_realigned.bam -q tumor_recal_data.table -q normal_recal_data.table \
    --algo Realigner -k $KNOWN_MILLS -k $KNOWN_INDELS tn_corealigned.bam
if [ "$?" -ne "0" ]; then   echo "Corealignment failed";   exit 1; fi

# ******************************************
# 7. Somatic Variant Calling
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i tn_corealigned.bam \
    --algo TNsnv --tumor_sample $TUMOR_SM --normal_sample $NORMAL_SM --dbsnp $KNOWN_DBSNP \
    --call_stats_out output-call.stats output-tnsnv.vcf.gz
if [ "$?" -ne "0" ]; then   echo "TNsnv failed";   exit 1; fi

$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i tn_corealigned.bam \
    --algo TNhaplotyper --tumor_sample $TUMOR_SM --normal_sample $NORMAL_SM \
    --dbsnp $KNOWN_DBSNP output-tnhaplotyper.vcf.gz
if [ "$?" -ne "0" ]; then   echo "TNhaplotyper failed";   exit 1; fi
