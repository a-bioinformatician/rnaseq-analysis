#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M ryanp_abo@dfci.harvard.edu
#$ -pe pvm 8
###################

# ------------------------------
# Script to run RNA-SeQC program
# ------------------------------

# SOFTWARE
rnaseqc='/ifs/rcgroups/ccgd/software/RNA-SeQC_v1.1.7.jar'
bwa='/ifs/rcgroups/ccgd/software/bwa-0.6.2/bwa' #7.5a/bwa'
picard='/ifs/rcgroups/ccgd/software/picard-tools-1.93'
samtools='/ifs/rcgroups/ccgd/rpa4/software/samtools-1.1/bin/samtools' #/usr/local/samtools-0.1.18-default/samtools'

BAM=$1
SAMPLE_NAME=$2
OUT_DIR=$3

REF_FASTA_BASE=/ifs/rcgroups/ccgd/reference/human/gencode/GRCh37-p13/genome/primary-assembly/GRCh37-p13.primary_assembly.genome #$REF_PATH/current/genome/primary-assembly/$REF_NAME
RNA_FASTA=/ifs/rcgroups/ccgd/reference/human/refseq/rRNA/human_all_rRNA.fasta # $REF_PATH/human_all_rRNA/human_all_rRNA.fasta
GTF=/ifs/rcgroups/ccgd/reference/human/gencode/GRCh37-p13/annotation/primary-assembly/gencode.v19.annotation.gtf # $REF_PATH/current/annotation/primary-assembly/$REF_NAME.gtf

mkdir -p $OUT_DIR/temp
TMP_DIR=$OUT_DIR/temp
cd $TMP_DIR

# Add read group information
echo java -Xmx2G -jar $picard/AddOrReplaceReadGroups.jar I=$BAM O=$TMP_DIR/alignments.rg.bam RGSM=$SAMPLE_NAME RGID=$SAMPLE_NAME RGLB=TruSeq RGPL=illumina RGPU=barcode VALIDATION_STRINGENCY=SILENT TMP_DIR=add_rg_tag_tmp SORT_ORDER=coordinate >> rg_tag.log 2>&1
java -Xmx2G -jar $picard/AddOrReplaceReadGroups.jar I=$BAM O=$TMP_DIR/alignments.rg.bam RGSM=$SAMPLE_NAME RGID=$SAMPLE_NAME RGLB=TruSeq RGPL=illumina RGPU=barcode VALIDATION_STRINGENCY=SILENT TMP_DIR=add_rg_tag_tmp SORT_ORDER=coordinate >> rg_tag.log 2>&1

# Make sure reference fasta is indexed
if [ ! -e $REF_FASTA_BASE.fa.fai ]; then
  $samtools faidx $REF_FASTA_BASE.fa
fi

if [ ! -e $REF_FASTA_BASE.dict ]; then
  # Create sequence dictionary
  echo java -jar $picard/CreateSequenceDictionary.jar REFERENCE=$REF_FASTA_BASE.fa OUTPUT=$REF_FASTA_BASE.dict
  java -jar $picard/CreateSequenceDictionary.jar REFERENCE=$REF_FASTA_BASE.fa OUTPUT=$REF_FASTA_BASE.dict
fi

echo java -jar $picard/ReorderSam.jar INPUT=$TMP_DIR/alignments.rg.bam OUTPUT=$TMP_DIR/alignments.rg.reorder.bam REFERENCE=$REF_FASTA_BASE.fa
java -jar $picard/ReorderSam.jar INPUT=$TMP_DIR/alignments.rg.bam OUTPUT=$TMP_DIR/alignments.rg.reorder.bam REFERENCE=$REF_FASTA_BASE.fa

# Mark duplicates
echo java -jar $picard/MarkDuplicates.jar INPUT=$TMP_DIR/alignments_rg.sorted2.bam OUTPUT=$TMP_DIR/alignments_rg.markeddups.bam METRICS_FILE=$TMP_DIR/duplicate_metrics ASSUME_SORTED=true 
java -jar $picard/MarkDuplicates.jar INPUT=$TMP_DIR/alignments.rg.reorder.bam OUTPUT=$TMP_DIR/alignments.rg.reorder.markeddups.bam METRICS_FILE=$TMP_DIR/duplicate_metrics ASSUME_SORTED=true 

# Index bam file
echo java -jar $picard/BuildBamIndex.jar INPUT=$TMP_DIR/alignments.rg.reorder.markeddups.bam
java -jar $picard/BuildBamIndex.jar INPUT=$TMP_DIR/alignments.rg.reorder.markeddups.bam

# RNA-SeQC
echo java -jar $rnaseqc -bwa $bwa -BWArRNA $RNA_FASTA -n 1000 -o $OUT_DIR -r $REF_FASTA_BASE.fa -s "$SAMPLE_NAME|$TMP_DIR/alignments.rg.reorder.markeddups.bam|$SAMPLE_NAME" -t $GTF
java -jar $rnaseqc -bwa $bwa -BWArRNA $RNA_FASTA -n 1000 -o $OUT_DIR -r $REF_FASTA_BASE.fa -s "$SAMPLE_NAME|$TMP_DIR/alignments.rg.reorder.markeddups.bam|$SAMPLE_NAME" -t $GTF

rm -r $TMP_DIR
