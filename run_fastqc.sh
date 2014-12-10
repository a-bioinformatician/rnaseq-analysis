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
SAMPLE_NAME=$1
READ_NUM=$2 # 1,2
FQ_DIR=$3
OUT_DIR=$4

mkdir -p $OUT_DIR
$FASTQC -t 4 -o $OUT_DIR $FQ_DIR/'read'$READ_NUM.fastq
