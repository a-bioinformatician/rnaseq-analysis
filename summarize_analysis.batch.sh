#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -pe pvm 4
###################

DIR=$1
SAMPLE_NAMES=$2

# SOFTWARE
summarize_script='/ifs/rcgroups/ccgd/rpa4/scripts/github_repos/rnaseq-analysis/summarize_analysis.sh'

cd $DIR

SAMPLE_NAMES_LST=$(echo $SAMPLE_NAMES | tr "," "\n")
echo 'Samples', $SAMPLE_NAMES_LST
for SAMPLE_NAME in $SAMPLE_NAMES_LST;
do
  echo $SAMPLE_NAME
  cd $SAMPLE_NAME
  qsub $summarize_script $DIR $DIR/$SAMPLE_NAME/summary $SAMPLE_NAME
  cd $DIR
done

