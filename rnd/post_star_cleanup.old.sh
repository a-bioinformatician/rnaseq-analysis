#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M ryanp_abo@dfci.harvard.edu
#$ -pe pvm 6
###################

# -----------------------------------------------
# Script to run STAR aligner
# ----------------------------------------------

samtools=/ifs/rcgroups/ccgd/rpa4/software/samtools-1.1/bin/samtools
OUT_DIR=$1 # Directory to place output files. This will be generated.
OUT_PREFIX=$2
REF_SRC=$3

mkdir -p $OUT_DIR

cd $OUT_DIR
if [ -e $OUT_PREFIX'_'$REF_SRC'_genome.bam' ] ; then
    rm $OUT_PREFIX'_'$REF_SRC'_genome.bam'
    rm $OUT_PREFIX'_'$REF_SRC'_trx.bam'
fi

# mv $OUT_PREFIX'_'$REF_SRC'_Aligned.sortedByCoord.out.bam' $OUT_PREFIX'.'$REF_SRC'_genome.bam'
# mv $OUT_PREFIX'_'$REF_SRC'_Aligned.toTranscriptome.out.bam' $OUT_PREFIX'.'$REF_SRC'_trx.bam'

if [ -e $OUT_PREFIX'_'$REF_SRC'_Aligned.sortedByCoord.out.bam' ] ; then
    mv $OUT_PREFIX'_'$REF_SRC'_Aligned.sortedByCoord.out.bam' $OUT_PREFIX'.'$REF_SRC'_genome.bam'
fi

if [ -e $OUT_PREFIX'_'$REF_SRC'_Aligned.toTranscriptome.out.bam' ] ; then
    mv $OUT_PREFIX'_'$REF_SRC'_Aligned.toTranscriptome.out.bam' $OUT_PREFIX'.'$REF_SRC'_trx.bam'
fi

$samtools sort $OUT_PREFIX'.'$REF_SRC'_genome.bam' $OUT_PREFIX'.'$REF_SRC'_genome.sorted'
echo $OUT_PREFIX'.'$REF_SRC'_genome.sorted.bam' $OUT_PREFIX'.'$REF_SRC'_genome.bam'
mv $OUT_PREFIX'.'$REF_SRC'_genome.sorted.bam' $OUT_PREFIX'.'$REF_SRC'_genome.bam'
