#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -pe pvm 2
###################

# -----------------------------------
# Script to run HTSeq for exon counts
# -----------------------------------

dexseq_count=/ifs/rcgroups/ccgd/rpa4/software/DEXSeq/inst/python_scripts/dexseq_count.py

ALIGNED_BAM=$1
GFF_FILE=$2
OUT_FILE=$3
LOG_FILE=$4

echo python $dexseq_count -p yes -s no -f bam -r pos $GFF_FILE $ALIGNED_BAM $OUT_FILE > $LOG_FILE
/usr/local/python-2.7.2/bin/python $dexseq_count -p yes -s no -f bam -r pos $GFF_FILE $ALIGNED_BAM $OUT_FILE
