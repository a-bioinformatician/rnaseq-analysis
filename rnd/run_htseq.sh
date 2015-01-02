#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -pe pvm 2
###################

# ------------------------------
# Script to run HTSeq
# ------------------------------

htseq_count=/ifs/rcgroups/ccgd/rpa4/bin/htseq-count

ALIGNED_BAM=$1
GTF_FILE=$2
OUT_FILE=$3
LOG_FILE=$4

echo python $htseq_count -f bam -r pos -m intersection-nonempty -s no $ALIGNED_BAM $GTF_FILE $OUT_FILE > $LOG_FILE
/usr/local/python-2.7.2/bin/python2.7 $htseq_count -f bam -r pos -m intersection-nonempty -s no $ALIGNED_BAM $GTF_FILE > $OUT_FILE
