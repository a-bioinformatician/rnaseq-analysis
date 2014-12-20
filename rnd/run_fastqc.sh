#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M ryanp_abo@dfci.harvard.edu
#$ -pe pvm 2
###################

# ------------------------------
# Script to run FastQC
# ------------------------------
# SOFTWARE
FASTQC='/ifs/rcgroups/ccgd/software/FastQC/fastqc'

# VARIABLES
FQ_DIR=$1
OUT_DIR=$2
LOG_FILE=$3

FQS=($(ls -d $FQ_DIR/*'.fastq.gz'))

mkdir -p $OUT_DIR
echo $FASTQC -t 4 -o $OUT_DIR $FQS >> $LOG_FILE
$FASTQC -t 4 -o $OUT_DIR $FQS
