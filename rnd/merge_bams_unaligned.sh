#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -pe pvm 2
###################

samtools=/ifs/rcgroups/ccgd/rpa4/software/samtools-1.1/bin/samtools
IN_BAM1=$1
IN_BAM2=$2
OUT_BAM=$3

cd $OUT_DIR

$samtools merge $OUT_BAM $IN_BAM1 $IN_BAM2
