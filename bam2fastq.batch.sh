#!/bin/sh

###############################################
# For an experiment or project with a set number
# of samples, each with an unaligned BAM file
# or multiple lanes of BAM files, convert each
# set to fastq files (paired-end).
# 
# Convert all bam files to fastq files.
###############################################

# SOFTWARE
BAM2FQ='/ifs/rcgroups/ccgd/rpa4/scripts/github_repos/rnaseq-analysis/bam2fastq.sh'

# SAMPLES
PROJECT_DIR=$1
SAMPLE_NAMES=$2

cd $PROJECT_DIR

SAMPLE_NAMES=(`echo $SAMPLE_NAMES | tr "," "\n"`)
echo 'Converting BAM to Fastq for' $SAMPLE_NAMES
for SAMPLE_NAME in "${SAMPLE_NAMES[@]}"
do
  mkdir -p $PWD/$SAMPLE_NAME/logs
  FASTQ_OUTDIR=$PWD/$SAMPLE_NAME/data
  FQ1=$FASTQ_OUTDIR/'read1.fastq'
  FQ2=$FASTQ_OUTDIR/'read2.fastq'

  SAMPLE_BAM_PATHS=`ls -d $FASTQ_OUTDIR/*bam | tr '\n' ','` # PATH TO SAMPLE BAM FILES.

  # GENERATE FASTQ DATA FROM ALIGNED SAMPLE BAM FILE
  if [ -e $FQ1 ] && [ -e $FQ2 ] ; then
    echo 'Fastq files previously generated'
  else 
    cd $SAMPLE_NAME/logs
    echo 'Generating fastq files from BAM file.'
    echo sh $BAM2FQ $SAMPLE_BAM_PATHS $FASTQ_OUTDIR
    qsub $BAM2FQ $SAMPLE_BAM_PATHS $FASTQ_OUTDIR
  fi
  cd $PROJECT_DIR
done
