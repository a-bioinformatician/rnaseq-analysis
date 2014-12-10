#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M ryanp_abo@dfci.harvard.edu
###################

# SAMPLES
SAMPLE_BAM_PATHS=$1 # PATHS TO SAMPLE BAM FILES. COMMA DELIMITED FOR MULTIPLE BAMS.
OUT_DIR=$2 # DIRECTORY TO STORE FASTQ FILES

bamutil='/ifs/rcgroups/ccgd/software/bamUtil-master/bin/bam'
BAMS=(`echo $SAMPLE_BAM_PATHS | tr "," "\n"`)

# GENERATE FASTQ DATA FROM ALIGNED SAMPLE BAM FILE
mkdir -p $OUT_DIR/tmp
cd $OUT_DIR
ITER=0
for BAM_F in "${BAMS[@]}"
do
#  echo /usr/local/java/bin/java -Dsamjdk.use_async_io=true -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx256m -Xms256m -jar /ifs/rcgroups/ccgd/software/picard-tools-1.92/SamToFastq.jar INPUT=$BAM_F FASTQ=$BAM_F.read1.fastq.gz SECOND_END_FASTQ=$BAM_F.read2.fastq.gz TMP_DIR=$OUT_DIR/tmp RE_REVERSE=true INCLUDE_NON_PF_READS=true CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 
#  /usr/local/java/bin/java -Dsamjdk.use_async_io=true -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx256m -Xms256m -jar /ifs/rcgroups/ccgd/software/picard-tools-1.92/SamToFastq.jar INPUT=$BAM_F FASTQ=$BAM_F.read1.fastq.gz SECOND_END_FASTQ=$BAM_F.read2.fastq.gz TMP_DIR=$OUT_DIR/tmp RE_REVERSE=true INCLUDE_NON_PF_READS=true CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2
#  gunzip $BAM_F.read1.fastq.gz
#  gunzip $BAM_F.read2.fastq.gz
#  if [ $ITER == 0 ]; then
#    mv $BAM_F.read1.fastq read1.fastq.tmp
#    mv $BAM_F.read2.fastq read2.fastq.tmp
#  else 
#    cat read1.fastq $BAM_F.read1.fastq > read1.fastq.tmp
#    cat read2.fastq $BAM_F.read2.fastq > read2.fastq.tmp
#    rm $BAM_F.read1.fastq
#    rm $BAM_F.read2.fastq
#  fi
#  mv read1.fastq.tmp read1.fastq
#  mv read2.fastq.tmp read2.fastq

  echo $bamutil bam2FastQ --in $BAM_F --outBase $BAM_F
  $bamutil bam2FastQ --in $BAM_F --outBase $BAM_F
  
  if [ $ITER == 0 ]; then
    mv $BAM_F'_1.fastq' read1.fastq.tmp
    mv $BAM_F'_2.fastq' read2.fastq.tmp
  else
    cat read1.fastq $BAM_F'_1.fastq' > read1.fastq.tmp
    cat read2.fastq $BAM_F'_2.fastq' > read2.fastq.tmp
    rm $BAM_F'_1.fastq'
    rm $BAM_F'_2.fastq'
  fi
  mv read1.fastq.tmp read1.fastq
  mv read2.fastq.tmp read2.fastq

  ITER=$(( $ITER + 1 ))
done
