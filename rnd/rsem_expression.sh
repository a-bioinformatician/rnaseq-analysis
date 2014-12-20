#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M ryanp_abo@dfci.harvard.edu
#$ -pe pvm 8
###################

###############################################
# RSEM - generates expression values 
###############################################

RSEM='/ifs/rcgroups/ccgd/rpa4/software/rsem-1.2.19'

SAMPLE_NAME=$1
SAMPLE_BAM=$2 # PATH TO TRANSCRIPTOME-ALIGNED BAM FILE.
RSEM_REF_PATH=$3 # PATH TO RSEM REFERENCE DIRECTORY WITH BASE_NAME
OUT=$4

mkdir -p $OUT
cd $OUT

echo $RSEM/rsem-calculate-expression -p 8 --bam --paired-end $SAMPLE_BAM $RSEM_REF_PATH $SAMPLE_NAME
$RSEM/rsem-calculate-expression -p 8 --bam --paired-end $SAMPLE_BAM $RSEM_REF_PATH $SAMPLE_NAME
