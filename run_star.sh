#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M ryanp_abo@dfci.harvard.edu
#$ -pe pvm 16
###################

# -----------------------------------------------
# Script to run STAR aligner
# ----------------------------------------------

# STAR='/ifs/rcgroups/ccgd/software/STAR_2.4.0h/bin/Linux_x86_64_static/STAR'
REF_BASE='/ifs/rcgroups/ccgd/reference'

FQ_DIR=$1 # Expecting paired-end fastq files named read1.fastq and read2.fastq
OUT_DIR=$2 # Directory to place output files. This will be generated.
SPECIES=$3 # human, mouse
REF_SRC=$4 # gencode, refseq
STAR_VER=$5 # 2.4.0h

STAR='/ifs/rcgroups/ccgd/software/STAR_'$STAR_VER'/bin/Linux_x86_64_static/STAR'

# Store the reference files in the CCGD reference directory for the current reference source.
REF_PATH=$(readlink -f $REF_BASE/$SPECIES/$REF_SRC/current)
STAR_REF=$REF_PATH/genome/star_$STAR_VER

FQ1=$FQ_DIR/'read1.fastq'
FQ2=$FQ_DIR/'read2.fastq'

if [ -e $FQ1 ] && [ -e $FQ2 ] ; then
  echo 'Fastq files previously generated'
else 
  echo 'Generate fastq files from unaligned BAM file before preceeding.'
  exit
fi

mkdir -p $OUT_DIR
cd $OUT_DIR


echo $STAR --genomeDir $STAR_REF --readFilesIn $FQ_DIR/'read1.fastq' $FQ_DIR/'read2.fastq' --runThreadN 16 --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --chimSegmentMin 20 --chimJunctionOverhangMin 20
$STAR --genomeDir $STAR_REF --readFilesIn $FQ_DIR/'read1.fastq' $FQ_DIR/'read2.fastq' --runThreadN 16 --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --chimSegmentMin 20 --chimJunctionOverhangMin 20

# echo $STAR --genomeDir $STAR_REF --readFilesIn $FQ_DIR/'read1.fastq' $FQ_DIR/'read2.fastq' --runThreadN 12 --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate --chimSegmentMin 20 --chimJunctionOverhangMin 20
# $STAR --genomeDir $STAR_REF --readFilesIn $FQ_DIR/'read1.fastq' $FQ_DIR/'read2.fastq' --runThreadN 12 --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate --chimSegmentMin 20 --chimJunctionOverhangMin 20
