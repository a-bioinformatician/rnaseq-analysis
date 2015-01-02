#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M ryanp_abo@dfci.harvard.edu
#$ -pe pvm 8
###################

# ------------------------------
# Script to prepare RSEM reference files. This can be run once and 
# reused.
# ------------------------------

RSEM='/ifs/rcgroups/ccgd/software/rsem-1.2.19'
BOWTIE='/ifs/rcgroups/ccgd/software/bowtie-1.0.0'
REF_BASE='/ifs/rcgroups/ccgd/reference'

SPECIES=$1 # human, mouse
REF_SRC=$2 # gencode, refseq
RSEM_VER=$3
LOG_FILE=$4

# Store the reference files in the CCGD reference directory for the current reference source.
REF_PATH=$(readlink -f $REF_BASE/$SPECIES/$REF_SRC/current)
GTF=$(readlink -f $REF_PATH/annotation/gtf)
REF_FASTA=$(readlink -f $REF_PATH/genome/fasta)
OUT=$REF_PATH/transcriptome/rsem_$RSEM_VER

mkdir -p $OUT
cd $OUT

OUT_BASENAME=$(basename $GTF '.gtf')

echo $RSEM/rsem-prepare-reference --bowtie-path $BOWTIE --gtf $GTF $REF_FASTA $OUT/$OUT_BASENAME > $LOG_FILE
$RSEM/rsem-prepare-reference --bowtie-path $BOWTIE --gtf $GTF $REF_FASTA $OUT/$OUT_BASENAME
