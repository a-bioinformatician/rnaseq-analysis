#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
###################

# ------------------------------
# Script to prepare RSEM reference
# files. This can be run once and 
# reused.
# ------------------------------

SPECIES=human

# SOFTWARE
BOWTIE='/ifs/rcgroups/ccgd/software/bowtie-1.0.0'
RSEM='/ifs/rcgroups/ccgd/software/rsem-1.2.19'

# PARAMS
REF_SRC=$1 # gencode
OUTDIR_BASE=$2 # <base directory for analysis>

REF_PATH='/ifs/rcgroups/ccgd/rpa4/ref2/'$SPECIES'/'$REF_SRC
REF_NAME=`cat $REF_PATH/current_version`
GTF=$REF_PATH'/current/annotation/gtf'
REF_FASTA=$REF_PATH'/current/genome/primary-assembly/split'
OUT=$OUTDIR_BASE/$REF_SRC/$REF_NAME

mkdir -p $OUTDIR_BASE/$REF_SRC
cd $OUTDIR_BASE/$REF_SRC
echo $RSEM/rsem-prepare-reference --bowtie-path $BOWTIE --no-polyA --gtf $GTF $REF_FASTA $OUT 
$RSEM/rsem-prepare-reference --bowtie-path $BOWTIE --no-polyA --gtf $GTF $REF_FASTA $OUT 
