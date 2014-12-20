#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M ryanp_abo@dfci.harvard.edu
#$ -pe pvm 15
###################

# -----------------------------------------------
# Script to run STAR aligner
# ----------------------------------------------

STAR='/ifs/rcgroups/ccgd/software/STAR_2.4.0f1/bin/Linux_x86_64_static/STAR'

FQ1=$1 # Expecting paired-end fastq files named read1.fastq and read2.fastq
FQ2=$2
OUT_DIR=$3 # Directory to place output files. This will be generated.
OUT_PREFIX=$4
STAR_REF=$5 # Full path to reference directory
LOG_FILE=$6

mkdir -p $OUT_DIR

echo $STAR --genomeDir $STAR_REF --readFilesCommand zcat --readFilesIn $FQ_DIR/'read1.fastq' $FQ_DIR/'read2.fastq' --runThreadN 15 --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --chimSegmentMin 20 --chimJunctionOverhangMin 20 --bamRemoveDuplicatesType UniqueIdential --outFileNamePrefix $OUT_DIR/$OUT_PREFIX'_' --genomeLoad LoadAndKeep>> $LOG_FILE
$STAR --genomeDir $STAR_REF --readFilesCommand zcat --readFilesIn $FQ1 $FQ2 --runThreadN 15 --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --chimSegmentMin 20 --chimJunctionOverhangMin 20 --bamRemoveDuplicatesType UniqueIdentical --outFileNamePrefix $OUT_DIR/$OUT_PREFIX'_' --genomeLoad LoadAndKeep

cd $OUT_DIR
ln -s $OUT_DIR/$OUT_PREFIX'_Aligned.sortedByCoord.out.bam' $OUT_PREFIX'_genome.bam'
ln -s $OUT_DIR/$OUT_PREFIX'_Aligned.toTranscriptome.out.bam' $OUT_PREFIX'_trx.bam'
