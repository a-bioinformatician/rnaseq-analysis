#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M ryanp_abo@dfci.harvard.edu
###################

# SAMPLES
BAMF=$1 # PATHS TO SAMPLE BAM FILES. COMMA DELIMITED FOR MULTIPLE BAMS.
SAMPLE_NAME=$2
BAM_OUTDIR=$3 # DIRECTORY TO STORE FASTQ FILES
METRICS_OUTDIR=$4
REMOVE_DUPS=$5

mkdir -p $BAM_OUTDIR/tmp
mkdir -p $METRICS_OUTDIR

ALIGNED_BASENAME=$(basename $BAMF .bam)

/usr/local/java/bin/java -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms2G  -Xmx4G -jar /ifs/rcgroups/ccgd/software/picard-tools-1.93/MarkDuplicates.jar INPUT=$BAMF OUTPUT=$BAM_OUTDIR/$ALIGNED_BASENAME'.dedupped.bam' METRICS_FILE=$METRICS_OUTDIR/$SAMPLE_NAME'_duplicateMetrics.txt' TMP_DIR=$BAM_OUTDIR/dedup_tmp MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=200 CREATE_INDEX=true CREATE_MD5_FILE=false READ_NAME_REGEX="[a-zA-Z0-9_]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*" COMPRESSION_LEVEL=1 VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=$REMOVE_DUPS
