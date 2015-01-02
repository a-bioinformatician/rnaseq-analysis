#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M ryanp_abo@dfci.harvard.edu
#$ -pe pvm 2
###################

FQ1=$1
FQ2=$2
N_DS=$3
OUT_DIR=$4
LOG_FILE=$5

echo $FQ1
echo $FQ2

FQ1_BASENAME=$(basename "$FQ1" .fastq.gz)
FQ2_BASENAME=$(basename "$FQ2" .fastq.gz)
echo $FQ1_BASENAME
echo $FQ2_BASENAME

export FQ1_DS=$OUT_DIR/$FQ1_BASENAME.fastq
export FQ2_DS=$OUT_DIR/$FQ2_BASENAME.fastq

mkdir -p $OUT_DIR/tmp
cd $OUT_DIR

gunzip -c $FQ1 > tmp/read1.fastq
gunzip -c $FQ2 > tmp/read2.fastq

FQ1=$PWD/tmp/read1.fastq
FQ2=$PWD/tmp/read2.fastq

paste $FQ1 $FQ2 | \
awk 'BEGIN{srand()}; {OFS="\t"; \
                      getline seqs; getline sep; gelint quals; \
                      print rand(),$0,seqs,sep,quals}' | \
sort -k1,1 | \
head -n $N_DS | \
awk '{OFS="\n"; \
      print $2,$4,$6,$8 >> ENVIRON["FQ1_DS"]; \
      print $3,$5,$7,$9 >> ENVIRON["FQ2_DS"]}'

#rm -r $OUT_DIR/tmp
