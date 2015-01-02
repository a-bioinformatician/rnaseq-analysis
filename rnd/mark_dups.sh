#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M ryanp_abo@dfci.harvard.edu
#$ -pe pvm 2
###################

markdups='/ifs/rcgroups/ccgd/software/picard-tools-1.93/MarkDuplicates.jar'
IN_BAM=$1 
METRICS_OUTFILE=$2
TMP_DIR=$3
LOG_FILE=$4

mkdir -p $TMP_DIR

mv $IN_BAM $IN_BAM'.tmp'
OUT_BAM=$IN_BAM

echo /usr/local/java/bin/java -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms2G  -Xmx4G -jar $markdups INPUT=$IN_BAM.tmp OUTPUT=$OUT_BAM METRICS_FILE=$METRICS_OUTFILE TMP_DIR=$TMP_DIR MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=200 CREATE_INDEX=true CREATE_MD5_FILE=false READ_NAME_REGEX="[a-zA-Z0-9_]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*" COMPRESSION_LEVEL=1 VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=false >> $LOG_FILE

/usr/local/java/bin/java -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms2G  -Xmx4G -jar $markdups INPUT=$IN_BAM.tmp OUTPUT=$OUT_BAM METRICS_FILE=$METRICS_OUTFILE TMP_DIR=$TMP_DIR MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=200 CREATE_INDEX=true CREATE_MD5_FILE=false READ_NAME_REGEX="[a-zA-Z0-9_]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*" COMPRESSION_LEVEL=1 VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=false

rm $IN_BAM.tmp
rm $OUT_BAM.bai
rm -r $TMP_DIR
