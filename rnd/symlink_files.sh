#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M ryanp_abo@dfci.harvard.edu
#$ -pe pvm 1
###################

# -----------------------------------------------
# Script to symlink files
# ----------------------------------------------

OUT_DIR=$1 # Directory to place output files. This will be generated.
BASENAME=$2
OUT_PREFIX=$3
LOG_FILE=$4

mkdir -p $OUT_DIR

cd $OUT_DIR
if [ -e $OUT_PREFIX'_genome.bam' ] ; then
    rm $OUT_PREFIX'_genome.bam'
fi

if [ -e $OUT_PREFIX'_trx.bam' ]; then
    rm $OUT_PREFIX'_trx.bam'
fi

ln -s $OUT_DIR/$BASENAME'_Aligned.sortedByCoord.out.bam' $OUT_PREFIX'_genome.bam'
ln -s $OUT_DIR/$BASENAME'_Aligned.toTranscriptome.out.bam' $OUT_PREFIX'_trx.bam'
