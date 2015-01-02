#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M ryanp_abo@dfci.harvard.edu
#$ -pe pvm 3
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

# Index genome bam
$samtools sort $OUT_PREFIX'.'$REF_SRC'_genome.bam' $OUT_PREFIX'.'$REF_SRC'_genome.sorted'
mv $OUT_PREFIX'.'$REF_SRC'_genome.sorted.bam' $OUT_PREFIX'.'$REF_SRC'_genome.bam'
$samtools index $OUT_PREFIX'.'$REF_SRC'_genome.bam'

# Convert chimeric sam file into sorted and indexed bam file
echo $samtools view -Shub $OUT_PREFIX'_'$REF_SRC'_Chimeric.out.sam' | $samtools sort - $OUT_PREFIX'_'$REF_SRC'_Chimeric.out' 
$samtools view -Shub $OUT_PREFIX'_'$REF_SRC'_Chimeric.out.sam' | $samtools sort - $OUT_PREFIX'_'$REF_SRC'_Chimeric.out' 
$samtools index $OUT_PREFIX'_'$REF_SRC'_Chimeric.out.bam' 
rm $OUT_PREFIX'_'$REF_SRC'_Chimeric.out.sam'
