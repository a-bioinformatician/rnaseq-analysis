#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M ryanp_abo@dfci.harvard.edu
#$ -pe pvm 12 
###################

# -----------------------------------------------
# Script to prepare STAR aligner reference files.
# This will generate a directory in the CCGD reference
# path.
# ----------------------------------------------

STAR='/ifs/rcgroups/ccgd/software/STAR_2.4.0g1/bin/Linux_x86_64_static/STAR'
REF_BASE='/ifs/rcgroups/ccgd/reference'

SPECIES=$1 # human, mouse
REF_SRC=$2 # gencode, refseq
STAR_VER=$3

# Store the reference files in the CCGD reference directory for the current reference source.
REF_PATH=$(readlink -f $REF_BASE/$SPECIES/$REF_SRC/current)
GTF=$(readlink -f $REF_PATH/annotation/gtf)
REF_FASTA=$(readlink -f $REF_PATH/genome/fasta)
OUT=$REF_PATH/genome/star_$STAR_VER

if [ -e $OUT ] ; then
  echo $OUT already exists! Generating STAR reference files.
fi

mkdir -p $OUT
cd $OUT
echo $STAR --runMode genomeGenerate --genomeDir $OUT --genomeFastaFiles $REF_FASTA --sjdbGTFfile $GTF --sjdbOverhang 99 --runThreadN 12
$STAR --runMode genomeGenerate --genomeDir $OUT --genomeFastaFiles $REF_FASTA --sjdbGTFfile $GTF --sjdbOverhang 99 --runThreadN 12
