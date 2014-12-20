#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M ryanp_abo@dfci.harvard.edu
#$ -pe pvm 4
###################

mergebams='/ifs/rcgroups/ccgd/software/picard-tools-1.93/MergeSamFiles.jar'
IN_BAMS=$1 
OUT_BAM=$2 
TMP_DIR=$3

echo /usr/local/java/bin/java -Xms2G  -Xmx16G -jar $mergebams $IN_BAMS OUTPUT=$OUT_BAM TMP_DIR=$TMP_DIR MAX_RECORDS_IN_RAM=5000000 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true
/usr/local/java/bin/java -Xms2G  -Xmx16G -jar $mergebams $IN_BAMS OUTPUT=$OUT_BAM TMP_DIR=$TMP_DIR MAX_RECORDS_IN_RAM=5000000 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true
