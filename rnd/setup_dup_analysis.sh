#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M ryanp_abo@dfci.harvard.edu
#$ -pe pvm 2
###################

markdups='/ifs/rcgroups/ccgd/software/picard-tools-1.93/MarkDuplicates.jar'
samtools=/ifs/rcgroups/ccgd/rpa4/software/samtools-1.1/bin/samtools

IN_BAM=$1 
DD_ALIGNED_BAM=$2
DD_READS_BAM=$3
REF_FASTA=$4
TMP_DIR=$5
LOG_FILE=$6

mkdir -p $TMP_DIR
METRICS_OUTFILE=$TMP_DIR/metrics.txt

echo /usr/local/java/bin/java -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms2G -Xmx4G -jar $markdups INPUT=$IN_BAM OUTPUT=$DD_ALIGNED_BAM METRICS_FILE=$METRICS_OUTFILE TMP_DIR=$TMP_DIR MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=200 CREATE_INDEX=true CREATE_MD5_FILE=false READ_NAME_REGEX="[a-zA-Z0-9_]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*" COMPRESSION_LEVEL=1 VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true >> $LOG_FILE

/usr/local/java/bin/java -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms2G -Xmx4G -jar $markdups INPUT=$IN_BAM OUTPUT=$DD_ALIGNED_BAM METRICS_FILE=$METRICS_OUTFILE TMP_DIR=$TMP_DIR MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=200 CREATE_INDEX=true CREATE_MD5_FILE=false READ_NAME_REGEX="[a-zA-Z0-9_]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*" COMPRESSION_LEVEL=1 VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

$samtools view $DD_ALIGNED_BAM | grep 'HI:i:1' | $samtools view -bT $REF_FASTA - > $DD_READS_BAM

rm $METRICS_OUTFILE
rm -r $TMP_DIR


