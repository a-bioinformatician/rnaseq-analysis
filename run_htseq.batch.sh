#!/bin/sh

###############################################
# Alignment QC software. 
###############################################

# SOFTWARE
SCRIPTS=/ifs/rcgroups/ccgd/rpa4/scripts/github_repos/rnaseq-analysis/

# SAMPLES
PROJECT_DIR=$1
SAMPLE_NAMES=$2
DEDUP=$3

cd $PROJECT_DIR

SAMPLE_NAMES=(`echo $SAMPLE_NAMES | tr "," "\n"`)
echo 'RNA-seq analysis for' $SAMPLE_NAMES
for SAMPLE_NAME in "${SAMPLE_NAMES[@]}"
do
  cd $SAMPLE_NAME
  mkdir -p htseq
  cd htseq

  BAMF=$PROJECT_DIR/$SAMPLE_NAME/star/Aligned.sortedByCoord.out.bam
  OUTF=htseq_refGene_all.txt
  if [ $DEDUP == 1 ]; then
    BAMF=$PROJECT_DIR/$SAMPLE_NAME/star/Aligned.sortedByCoord.out.dedupped.bam
    OUTF=htseq_refGene_all.dedupped.txt
  fi

  qsub -q ccgd_testq,all.q $SCRIPTS/run_htseq.sh $BAMF ~/ref2/human/refseq/current/annotation/hg19_refGene.gtf $OUTF

  cd $PROJECT_DIR
done
