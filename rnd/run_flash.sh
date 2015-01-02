#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -pe pvm 8
###################

# -----------------------------------
# Script to run FLASH
# -----------------------------------

flash=/ifs/rcgroups/ccgd/rpa4/bin/flash

FQ1=$1
FQ2=$2
OUT_BASENAME=$3
OUT_DIR=$4
LOG_FILE=$5

cd $OUT_DIR

echo $flash -t 8 -z -d $OUT_DIR -o $OUT_BASENAME -r 101 -M 80 $FQ1 $FQ2 > $LOG_FILE
$flash -t 8 -z -d $OUT_DIR -o $OUT_BASENAME -r 101 -M 80 $FQ1 $FQ2 &> flash.log 
