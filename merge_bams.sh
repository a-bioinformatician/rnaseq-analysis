#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -pe pvm 4
###################

OUT_DIR=$1
SAMPLE_NAME=$2
BAM_F1=$3
BAM_F2=$4

cd $OUT_DIR

samtools merge $SAMPLE_NAME.bam $BAM_F1 $BAM_F2
