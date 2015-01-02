#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M ryanp_abo@dfci.harvard.edu
#$ -pe pvm 4
###################

# ------------------------------
# Script to run FastQC
# ------------------------------
# SOFTWARE
FASTQC='/ifs/rcgroups/ccgd/software/FastQC/fastqc'

# VARIABLES
FQ_FILE=$1
OUT_DIR=$2
LOG_FILE=$3

mkdir -p $OUT_DIR
echo $FASTQC -t 4 -o $OUT_DIR $FQ_FILE >> $LOG_FILE
$FASTQC -t 4 -o $OUT_DIR $FQ_FILE
