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

IN_BAM=$1
OUT_BASENAME=$2
OUT_DIR=$3
REF_PATH=$4 # PATH TO RSEM REFERENCE DIRECTORY WITH BASE_NAME
LOG_FILE=$5

mkdir -p $OUT_DIR
cd $OUT_DIR

echo $RSEM/rsem-calculate-expression -p 8 --bam --paired-end $IN_BAM $REF_PATH $OUT_BASENAME
$RSEM/rsem-calculate-expression -p 8 --bam --paired-end $IN_BAM $REF_PATH $OUT_BASENAME
# $RSEM/rsem-calculate-expression -p 8 --calc-ci --bam --paired-end $IN_BAM $REF_PATH $OUT_BASENAME
