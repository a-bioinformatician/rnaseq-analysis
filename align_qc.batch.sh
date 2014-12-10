#!/bin/sh

###############################################
# Alignment QC software. 
###############################################

# SOFTWARE
SCRIPTS=/ifs/rcgroups/ccgd/rpa4/scripts/github_repos/rnaseq-analysis/
RSEQC=$SCRIPTS/rseq_qc.sh
RNASEQ_QC=$SCRIPTS/run_rnaseq_qc.sh

# SAMPLES
PROJECT_DIR=$1
SAMPLE_NAMES=$2

cd $PROJECT_DIR

SAMPLE_NAMES=(`echo $SAMPLE_NAMES | tr "," "\n"`)
echo 'RNA-seq analysis for' $SAMPLE_NAMES
for SAMPLE_NAME in "${SAMPLE_NAMES[@]}"
do
  LOG_DIR=$PWD/$SAMPLE_NAME/logs
  mkdir -p $LOG_DIR
  
  GENOME_ALIGN_BAM=$PROJECT_DIR/$SAMPLE_NAME/star/Aligned.out.sorted.bam

  FASTQ_OUTDIR=$PWD/$SAMPLE_NAME/data
  FQ1=$FASTQ_OUTDIR/'read1.fastq'
  FQ2=$FASTQ_OUTDIR/'read2.fastq'

  cd $LOG_DIR

  # RUN RSEQC 
  RSEQC_OUTDIR=$PROJECT_DIR/$SAMPLE_NAME/qc
  echo $RSEQC $SAMPLE_NAME $GENOME_ALIGN_BAM $RSEQC_OUTDIR/rseqc
  qsub -q ccgd.q,all.q $RSEQC $SAMPLE_NAME $GENOME_ALIGN_BAM $RSEQC_OUTDIR/rseqc

  # RUN RNASEQ QC
  RNASEQ_QC_OUTDIR=$PROJECT_DIR/$SAMPLE_NAME/qc
  echo $RNASEQ_QC $GENOME_ALIGN_BAM $SAMPLE_NAME $RNASEQ_QC_OUTDIR/rnaseqc
  qsub -q ccgd.q,all.q $RNASEQ_QC $GENOME_ALIGN_BAM $SAMPLE_NAME $RNASEQ_QC_OUTDIR/rnaseqc

  cd $PROJECT_DIR
done
