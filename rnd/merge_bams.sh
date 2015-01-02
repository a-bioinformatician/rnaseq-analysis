#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M ryanp_abo@dfci.harvard.edu
#$ -pe pvm 2
###################

markdups='/ifs/rcgroups/ccgd/software/picard-tools-1.93/MarkDuplicates.jar'
IN_BAM1=$1
IN_BAM2=$2
OUT_BAM=$3 
METRICS_OUTFILE=$4
TMP_DIR=$5
LOG_FILE=$6

mkdir -p $TMP_DIR

echo /usr/local/java/bin/java -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms2G  -Xmx4G -jar $markdups INPUT=$IN_BAM1 INPUT=$IN_BAM2 OUTPUT=$OUT_BAM METRICS_FILE=$METRICS_OUTFILE TMP_DIR=$TMP_DIR MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=200 CREATE_INDEX=true CREATE_MD5_FILE=false READ_NAME_REGEX="[a-zA-Z0-9_]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*" COMPRESSION_LEVEL=1 VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=false >> $LOG_FILE

/usr/local/java/bin/java -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms2G  -Xmx4G -jar $markdups INPUT=$IN_BAM1 INPUT=$IN_BAM2 OUTPUT=$OUT_BAM METRICS_FILE=$METRICS_OUTFILE TMP_DIR=$TMP_DIR MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=200 CREATE_INDEX=true CREATE_MD5_FILE=false READ_NAME_REGEX="[a-zA-Z0-9_]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*" COMPRESSION_LEVEL=1 VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=false

rm -r $TMP_DIR
