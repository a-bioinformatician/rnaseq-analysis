#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M ryanp_abo@dfci.harvard.edu
#$ -pe pvm 8
###################

###############################################
# RSEM EXPRESSION - generates expression values 
# using RSEM. 
#
# Ryan Abo
# July 23, 2014
#
# Requires RSEM reference prior to running this
#
###############################################

# SOFTWARE
RSEM='/ifs/rcgroups/ccgd/software/rsem-1.2.19'
BAM2FQ='/ifs/rcgroups/ccgd/rpa4/scripts/github_repos/rnaseq-analysis/bam2fastq.sh'
BOWTIE='/ifs/rcgroups/ccgd/software/bowtie-1.0.0'

# SAMPLES
SAMPLE_NAME=$1
SAMPLE_BAM_PATHS=$2 # PATH TO SAMPLE BAM FILES.
FASTQ_OUTDIR=$3 # DIRECTORY TO STORE FASTQ FILES
RSEM_REF_PATH=$4 # PATH TO THE RSEM REFERENCE FILES
OUT_DIR=$5

FQ1=$FASTQ_OUTDIR/'read1.fastq'
FQ2=$FASTQ_OUTDIR/'read2.fastq'

# GENERATE FASTQ DATA FROM ALIGNED SAMPLE BAM FILE
if [ -e $FQ1 ] && [ -e $FQ2 ] ; then
  echo 'Fastq files previously generated'
else 
  echo 'Generating fastq files from BAM file.'
  echo sh $BAM2FQ $SAMPLE_BAM_PATHS $FASTQ_OUTDIR
  sh $BAM2FQ $SAMPLE_BAM_PATHS $FASTQ_OUTDIR
fi

mkdir -p $OUT_DIR
cd $OUT_DIR
echo $RSEM/rsem-calculate-expression -p 8 --paired-end --bowtie-path $BOWTIE --bowtie-chunkmbs 1024 $FQ1 $FQ2 $RSEM_REF_PATH $SAMPLE_NAME
$RSEM/rsem-calculate-expression -p 8 --paired-end --bowtie-path $BOWTIE --bowtie-chunkmbs 1024 $FQ1 $FQ2 $RSEM_REF_PATH $SAMPLE_NAME

# Removing bam files corresponding to bowtie alignment to transcriptome.
# rm $OUT_DIR/*bam
