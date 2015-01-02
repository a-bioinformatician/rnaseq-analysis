#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -pe pvm 4
###################

# ------------------------------
# Script to run BWA
# ------------------------------

FQ1=$1
FQ2=$2
OUT_DIR=$3
LOG_FILE=$4

# SOFTWARE
bwa='/ifs/rcgroups/ccgd/software/bwa-0.7.5a'
samtools='/ifs/rcgroups/ccgd/rpa4/software/samtools-1.1/bin/samtools'

# PARAMS
mkdir -p $OUT_DIR

GENOME='/ifs/rcgroups/ccgd/reference/human/refseq/rRNA/human_all_rRNA.fasta'

FQ_BASENAME=$(basename $FQ1 '_1.fastq.gz')

# if [ -e $OUT_DIR/$FQ_BASENAME'.ribo_aligned.bam' ] ; then
#    echo 'Ribo aligned BAM previously generated for sample fastqs' $FQ_BASENAME >> $LOG_FILE
# else
#    $bwa/bwa mem -t 4 $GENOME $FQ1 $FQ2 | $samtools view -Shub - | $samtools sort - $OUT_DIR/$FQ_BASENAME'.ribo_aligned'
#    $samtools index $OUT_DIR/$FQ_BASENAME'.ribo_aligned.bam'
# fi

# $samtools flagstat $OUT_DIR/$FQ_BASENAME'.ribo_aligned.bam' > $OUT_DIR/$FQ_BASENAME'_ribo_aligned.flagstat.txt'
$samtools idxstats $OUT_DIR/$FQ_BASENAME'.ribo_aligned.bam' > $OUT_DIR/$FQ_BASENAME'_ribo_aligned.idxstats.txt'

# Do not remove yet
## rm $OUT_DIR/*'.bam'*
