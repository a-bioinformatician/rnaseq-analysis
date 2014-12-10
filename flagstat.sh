#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -pe pvm 4
###################

# ------------------------------
# Script to run samtools flagstat
# ------------------------------

BAM=$1
OUT_FILE=$2

echo samtools flagstat $BAM '>' $OUT_FILE
samtools flagstat $BAM > $OUT_FILE
