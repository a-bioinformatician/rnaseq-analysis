#!/bin/sh

###############################################
# For an experiment or project with a set number
# of samples, each with an unaligned BAM file
# or multiple lanes of BAM files, convert each
# set to fastq files (paired-end).
# 
# Convert all bam files to fastq files.
###############################################

SCRIPTS_DIR='/ifs/rcgroups/ccgd/rpa4/scripts/github_repos/rnaseq-analysis/rnd'
ALIGN_RRNA=$SCRIPTS_DIR/'align_rrna.sh'

PROJECT_DIR=$1
DATA_TAG=$2 # Name of directory containing fastq files to align.
SAMPLE_IDS=$3

cd $PROJECT_DIR

SAMPLE_ID_ARR=(`echo $SAMPLE_IDS | tr "," "\n"`)
echo 'Aligning rRNA for samples:' $SAMPLE_ID_ARR
for SAMPLE_ID in "${SAMPLE_ID_ARR[@]}"
do
    LOG_DIR=$PROJECT_DIR/$SAMPLE_ID/logs
    mkdir -p $LOG_DIR

    LANE_FILES=($(ls -d $PROJECT_DIR/$SAMPLE_ID/data/$DATA_TAG/$SAMPLE_ID*'.bam'))
    for LANE_FILE in "${LANE_FILES[@]}"
    do
        LANE_BASENAME=$(basename $LANE_FILE '.bam')
        echo $LANE_BASENAME
        FQ1=$PROJECT_DIR/$SAMPLE_ID/data/$DATA_TAG/$LANE_BASENAME'_1.fastq.gz'
        FQ2=$PROJECT_DIR/$SAMPLE_ID/data/$DATA_TAG/$LANE_BASENAME'_2.fastq.gz'
        OUTDIR=$PROJECT_DIR/$SAMPLE_ID/align/ribo
        cd $SAMPLE_ID/logs
        echo $0 'Aligning reads to rRNA for sample' $SAMPLE_ID 'fastq files in' $FQ_DIR >> $LOG_DIR/analysis.log
        echo 'qsub' $ALIGN_RRNA $FQ1 $FQ2 $OUTDIR $LOG_DIR/analysis.log >> $LOG_DIR/analysis.log
        qsub $ALIGN_RRNA $FQ1 $FQ2 $OUTDIR $LOG_DIR/analysis.log
    done 
    cd $PROJECT_DIR
done
