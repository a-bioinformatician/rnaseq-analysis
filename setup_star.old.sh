#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M ryanp_abo@dfci.harvard.edu
#$ -pe pvm 8 
###################

# ------------------------------
# Script to prepare STAR aligner
# ------------------------------

SPECIES=human

# SOFTWARE
STAR='/ifs/rcgroups/ccgd/software/STAR_2.3.0e/STAR'

# PATHS
BASE='/ifs/rcgroups/ccgd/rpa4/analysis/rnaseq/130701'

# PARAMS
REF_SRC=$1
REF_PATH='/ifs/rcgroups/ccgd/rpa4/ref2/'$SPECIES'/'$REF_SRC
REF_NAME=`cat $REF_PATH/current_version`
GTF=$REF_PATH'/current/annotation/gtf'
REF_FASTA=$REF_PATH'/current/genome/fasta'
OUT=$BASE/$REF_SRC/star/ref

mkdir -p $OUT
cd $OUT
echo $STAR --runMode genomeGenerate --genomeDir $OUT --genomeFastaFiles $REF_FASTA --sjdbOverhang 100 --runThreadN 8
$STAR --runMode genomeGenerate --genomeDir $OUT --genomeFastaFiles $REF_FASTA --sjdbOverhang 100 --runThreadN 8
