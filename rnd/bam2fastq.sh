#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M ryanp_abo@dfci.harvard.edu
#$ -pe pvm 2
###################

SAMPLE_BAM_PATH=$1 # Path to sample bam file
OUTBASE=$2 # Directory and basename of the output fastq files.
LOG_FILE=$3

echo $LOG_FILE
bamutil='/ifs/rcgroups/ccgd/software/bamUtil-1.0.7/bin/bam'

if [ -e $OUTBASE'_1.fastq' ] && [ -e $OUTBASE'_2.fastq' ] ; then
    echo $0 'Fastq files previously generated for sample' $SAMPLE_ID 'for data tag' $DATA_TAG >> $LOG_FILE
elif [ -e $OUTBASE'_1.fastq.gz' ] && [ -e $OUTBASE'_2.fastq.gz' ] ; then
    echo $0 'Fastq files previously generated for sample' $SAMPLE_ID 'for data tag' $DATA_TAG >> $LOG_FILE
else 
    echo $bamutil bam2FastQ --in $SAMPLE_BAM_PATH --outBase $SAMPLE_BAM_PATH >> $LOG_FILE
    $bamutil bam2FastQ --in $SAMPLE_BAM_PATH --outBase $OUTBASE
fi

if [ -e $OUTBASE'_1.fastq' ] ; then 
    gzip $OUTBASE'_1.fastq'
fi

if [ -e $OUTBASE'_2.fastq' ] ; then
    gzip $OUTBASE'_2.fastq'
fi

# Old code using picard tools. It was buggy.
#  echo /usr/local/java/bin/java -Dsamjdk.use_async_io=true -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx256m -Xms256m -jar /ifs/rcgroups/ccgd/software/picard-tools-1.92/SamToFastq.jar INPUT=$BAM_F FASTQ=$BAM_F.read1.fastq.gz SECOND_END_FASTQ=$BAM_F.read2.fastq.gz TMP_DIR=$OUT_DIR/tmp RE_REVERSE=true INCLUDE_NON_PF_READS=true CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 
#  /usr/local/java/bin/java -Dsamjdk.use_async_io=true -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx256m -Xms256m -jar /ifs/rcgroups/ccgd/software/picard-tools-1.92/SamToFastq.jar INPUT=$BAM_F FASTQ=$BAM_F.read1.fastq.gz SECOND_END_FASTQ=$BAM_F.read2.fastq.gz TMP_DIR=$OUT_DIR/tmp RE_REVERSE=true INCLUDE_NON_PF_READS=true CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2
#  gunzip $BAM_F.read1.fastq.gz
#  gunzip $BAM_F.read2.fastq.gz
#  if [ $ITER == 0 ]; then
#    mv $BAM_F.read1.fastq read1.fastq.tmp
#    mv $BAM_F.read2.fastq read2.fastq.tmp
#  else 
#    cat read1.fastq $BAM_F.read1.fastq > read1.fastq.tmp
#    cat read2.fastq $BAM_F.read2.fastq > read2.fastq.tmp
#    rm $BAM_F.read1.fastq
#    rm $BAM_F.read2.fastq
#  fi
#  mv read1.fastq.tmp read1.fastq
#  mv read2.fastq.tmp read2.fastq


