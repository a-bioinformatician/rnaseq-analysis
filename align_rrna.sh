#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -pe pvm 4
###################

# ------------------------------
# Script to run BWA
# ------------------------------

FQ_DIR=$1
OUT_DIR=$2

# SOFTWARE
bwa='/ifs/rcgroups/ccgd/software/bwa-0.7.5a'
samtools='/usr/local/samtools-0.1.18-default'

# PARAMS
mkdir -p $OUT_DIR
GENOME='/ifs/rcgroups/ccgd/rpa4/ref2/human/rRNA/hgRibo'
OUT=$OUT_DIR/'hgRibo_align'
FQ1=$FQ_DIR/'read1.fastq'
FQ2=$FQ_DIR/'read2.fastq'

$bwa/bwa mem -t 4 $GENOME $FQ1 $FQ2 > $OUT_DIR/'hgRibo_align.sam'

$samtools/samtools view -Sb $OUT_DIR/'hgRibo_align.sam' > $OUT_DIR/'hgRibo_align.bam'
$samtools/samtools sort $OUT_DIR/'hgRibo_align.bam' $OUT_DIR/'hgRibo_align.sorted'
$samtools/samtools index $OUT_DIR/'hgRibo_align.sorted.bam'
$samtools/samtools flagstat $OUT_DIR/'hgRibo_align.sorted.bam' > $OUT_DIR/hgRibo_align_flagstat.txt
rm $OUT_DIR/'hgRibo_align.sam'
rm $OUT_DIR/'hgRibo_align.bam'
