#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M ryanp_abo@dfci.harvard.edu
#$ -pe pvm 2
###################

# -----------------------------------------------
# Script to sort bam by name
# ----------------------------------------------

samtools=/ifs/rcgroups/ccgd/rpa4/software/samtools-1.1/bin/samtools
IN_BAM=$1 # Directory to place output files. This will be generated.
LOG_FILE=$2

echo $samtools sort -n $IN_BAM $IN_BAM'.sortedname' >> $LOG_FILE
$samtools sort -n $IN_BAM $IN_BAM'.sortedname'
mv $IN_BAM'.sortedname.bam' $IN_BAM 
