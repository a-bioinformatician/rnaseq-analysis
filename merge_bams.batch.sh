#!/bin/sh

SCRIPTS=/ifs/rcgroups/ccgd/rpa4/scripts/github_repos/rnaseq-analysis/
MERGE=$SCRIPTS/merge_bams.sh

PROJECT_DIR=$1
SAMPLE_NAMES=$2

cd $PROJECT_DIR

SAMPLE_NAMES=(`echo $SAMPLE_NAMES | tr "," "\n"`)
for SAMPLE_NAME in "${SAMPLE_NAMES[@]}"
do

  BAM_F1=$PWD/$SAMPLE_NAME/data/$SAMPLE_NAME-lane1.bam
  BAM_F2=$PWD/$SAMPLE_NAME/data/$SAMPLE_NAME-lane2.bam

  echo $BAM_F1

  mkdir -p $SAMPLE_NAME/data
  cd $SAMPLE_NAME/data

  qsub $MERGE $PWD $SAMPLE_NAME $BAM_F1 $BAM_F2

  cd $PROJECT_DIR
done
