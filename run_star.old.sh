#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M ryanp_abo@dfci.harvard.edu
#$ -pe pvm 8 
###################

# ------------------------------
# Script to run STAR aligner
# ------------------------------

SPECIES=human

# SOFTWARE
STAR='/ifs/rcgroups/ccgd/software/STAR_2.3.1n/STAR'
SCRIPTS='/ifs/rcgroups/ccgd/rpa4/scripts'
TIME='/ifs/rcgroups/ccgd/rpa4/software/time-1.7/bin/time'

FQ_DIR=$1
OUT_DIR=$2

REF_PATH='/ifs/rcgroups/ccgd/rpa4/ref2/human/gencode'
GTF=$REF_PATH'/current/annotation/gtf'
STAR_REF=/ifs/rcgroups/ccgd/rpa4/analysis/rnaseq/gencode_ref/star/ref

mkdir -p $OUT_DIR
cd $OUT_DIR
echo $STAR --genomeDir $STAR_REF --readFilesIn $FQ_DIR/'read1.fastq' $FQ_DIR/'read2.fastq' --runThreadN 8 --sjdbGTFfile $GTF --chimSegmentMin 20 --chimJunctionOverhangMin 20
$TIME -v $STAR --genomeDir $STAR_REF --readFilesIn $FQ_DIR/'read1.fastq' $FQ_DIR/'read2.fastq' --runThreadN 8 --sjdbGTFfile $GTF --chimSegmentMin 20 --chimJunctionOverhangMin 20

sh $SCRIPTS/ccgd-utils/sort_index_bam.sh $OUT_DIR/Aligned.out.sam 1
sh $SCRIPTS/ccgd-utils/sort_index_bam.sh $OUT_DIR/Chimeric.out.sam 1
