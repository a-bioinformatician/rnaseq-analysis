#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -pe pvm 8
###################

####################################################
# PRADA pipeline - generates QC metrics and fusion
# calls.
#
# Ryan Abo
# July 11, 2014
#
####################################################

# SOFTWARE
PICARD='/ifs/rcgroups/ccgd/software/picard-tools-1.93'
BWA='/ifs/rcgroups/ccgd/rpa4/software/pyPRADA_1.2/tools/bwa-0.5.7-mh/bwa'
GATK='/ifs/rcgroups/ccgd/rpa4/software/pyPRADA_1.2/tools/GATK'
QC='/ifs/rcgroups/ccgd/software/RNA-SeQC_v1.1.7.jar'
PRADA_FUSION='/ifs/rcgroups/ccgd/rpa4/software/pyPRADA_1.2/prada-fusion'
SAMTOOLS='/ifs/rcgroups/ccgd/rpa4/software/samtools-1.1/bin/samtools'

#REFERENCE FILES
TRX_G_FA=/ifs/rcgroups/ccgd/rpa4/software/pyPRADA_1.2/PRADA-ref-hg19/Ensembl64.transcriptome.plus.genome.fasta
TRX_G_MAP=/ifs/rcgroups/ccgd/rpa4/software/pyPRADA_1.2/PRADA-ref-hg19/Ensembl64.transcriptome.plus.genome.map
HG19_FA=/ifs/rcgroups/ccgd/rpa4/software/pyPRADA_1.2/PRADA-ref-hg19/Homo_sapiens_assembly19.fasta
DBSNP_VCF=/ifs/rcgroups/ccgd/rpa4/software/pyPRADA_1.2/PRADA-ref-hg19/dbsnp_135.b37.vcf
HG19_GTF=/ifs/rcgroups/ccgd/rpa4/software/pyPRADA_1.2/PRADA-ref-hg19/Homo_sapiens.GRCh37.64.gtf
PRADA_CONF_F=/ifs/rcgroups/ccgd/rpa4/software/pyPRADA_1.2/conf.txt

# SAMPLE INFORMATION
SAMPLE_NAME=$1
FQ1=$2
FQ2=$3
OUT_DIR=$4

TMP_DIR=$OUT_DIR/temp
mkdir -p $TMP_DIR
cd $OUT_DIR
READS=( $FQ1 $FQ2 )
for FQ in "${READS[@]}"
do
  
  FQ_BASENAME=$(basename "$FQ" .fastq.gz)
  gunzip -c $FQ > $TMP_DIR/$FQ_BASENAME'.fastq'
  NEW_FQ=$TMP_DIR/$FQ_BASENAME'.fastq'

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # REALIGN THE READS TO THE TRANSCRIPTOME AND GENOME.
  echo $BWA aln -t 8 $TRX_G_FA $NEW_FQ > $TMP_DIR/$FQ_BASENAME'.sai'
  $BWA aln -t 8 $TRX_G_FA $NEW_FQ > $TMP_DIR/$FQ_BASENAME'.sai'

  echo $BWA samse -s -n 100 $TRX_G_FA $TMP_DIR/$FQ_BASENAME'.sai' $NEW_FQ > $TMP_DIR/$FQ_BASENAME'.sam'
  $BWA samse -s -n 100 $TRX_G_FA $TMP_DIR/$FQ_BASENAME'.sai' $NEW_FQ > $TMP_DIR/$FQ_BASENAME'.sam'
#  rm $TMP_DIR/$FQ_BASENAME'.sai'

  echo $SAMTOOLS view -bS -o $FQ_BASENAME'.bam' $TMP_DIR/$FQ_BASENAME'.sam'
  $SAMTOOLS view -bS -o $FQ_BASENAME'.bam' $TMP_DIR/$FQ_BASENAME'.sam'
#  rm $TMP_DIR/$FQ_BASENAME'.sam'

#  echo $SAMTOOLS sort -n -m 1000000000 $FQ_BASENAME'.bam' $FQ_BASENAME'.sorted'
#  $SAMTOOLS sort -n -m 1000000000 $FQ_BASENAME'.bam' $FQ_BASENAME'.sorted'
#  rm $FQ_BASENAME'.bam'
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # REMAP WITH GATK
  echo java -Djava.io.tmpdir=tmp/ -cp $GATK/RemapAlignments.jar -Xmx8g org.broadinstitute.cga.tools.gatk.rna.RemapAlignments M=$TRX_G_MAP IN=$FQ_BASENAME'.bam' OUT=$FQ_BASENAME'.remapped.bam' R=$HG19_FA REDUCE=TRUE
  java -Djava.io.tmpdir=tmp/ -cp $GATK/RemapAlignments.jar -Xmx8g org.broadinstitute.cga.tools.gatk.rna.RemapAlignments M=$TRX_G_MAP IN=$FQ_BASENAME'.bam' OUT=$FQ_BASENAME'.remapped.bam' R=$HG19_FA REDUCE=TRUE
#  rm -f $FQ_BASENAME'.sorted.bam'

  echo $SAMTOOLS sort -n -m 1000000000 $FQ_BASENAME'.remapped.bam' $FQ_BASENAME'.remapped.sorted'
  $SAMTOOLS sort -n -m 1000000000 $FQ_BASENAME'.remapped.bam' $FQ_BASENAME'.remapped.sorted'
#  rm -f $FQ_BASENAME'.remapped.bam'
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REJOIN PAIRED READS INTO SINGLE BAM FILE
FQ1_BASENAME=$(basename "$FQ1" .fastq.gz)
FQ2_BASENAME=$(basename "$FQ2" .fastq.gz)
echo java -Djava.io.tmpdir=tmp/ -Xmx8g -jar $GATK/PairMaker.jar IN1=$FQ1_BASENAME'.remapped.sorted.bam' IN2=$FQ2_BASENAME'.remapped.sorted.bam' OUTPUT=$SAMPLE_NAME.paired.bam TMP_DIR=tmp/
java -Djava.io.tmpdir=tmp/ -Xmx8g -jar $GATK/PairMaker.jar IN1=$FQ1_BASENAME'.remapped.sorted.bam' IN2=$FQ2_BASENAME'.remapped.sorted.bam' OUTPUT=$SAMPLE_NAME.paired.bam TMP_DIR=tmp/
#rm -f $FQ1_BASENAME'.remapped.sorted.bam'
#rm -f $FQ2_BASENAME'.remapped.sorted.bam'

echo $SAMTOOLS sort -m 1000000000 $SAMPLE_NAME.paired.bam $SAMPLE_NAME.paired.sorted
$SAMTOOLS sort -m 1000000000 $SAMPLE_NAME.paired.bam $SAMPLE_NAME.paired.sorted
#rm -f $SAMPLE_NAME.paired.bam
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
echo java -Xmx8g -jar $PICARD/AddOrReplaceReadGroups.jar I=$SAMPLE_NAME.paired.sorted.bam O=$SAMPLE_NAME.withRG.paired.sorted.bam RGLB=$SAMPLE_NAME RGPL=Illumina RGPU=$SAMPLE_NAME RGSM=$SAMPLE_NAME
java -Xmx8g -jar $PICARD/AddOrReplaceReadGroups.jar I=$SAMPLE_NAME.paired.sorted.bam O=$SAMPLE_NAME.withRG.paired.sorted.bam RGLB=$SAMPLE_NAME RGPL=Illumina RGPU=$SAMPLE_NAME RGSM=$SAMPLE_NAME
#rm -f $SAMPLE_NAME.paired.sorted.bam

echo $SAMTOOLS index $SAMPLE_NAME.withRG.paired.sorted.bam
$SAMTOOLS index $SAMPLE_NAME.withRG.paired.sorted.bam
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar -l INFO -R $HG19_FA --default_platform illumina --knownSites $DBSNP_VCF -I $SAMPLE_NAME.withRG.paired.sorted.bam --downsample_to_coverage 10000 -T CountCovariates -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -nt 8 -recalFile $SAMPLE_NAME.orig.csv

java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar -l INFO -R $HG19_FA --default_platform illumina -I $SAMPLE_NAME.withRG.paired.sorted.bam -T TableRecalibration --out $SAMPLE_NAME.withRG.GATKRecalibrated.bam -recalFile $SAMPLE_NAME.orig.csv
#rm -f $SAMPLE_NAME.withRG.paired.sorted.bam
#rm -f $SAMPLE_NAME.withRG.paired.sorted.bam.bai
#rm -f $SAMPLE_NAME.orig.csv

java -Xmx8g -jar $PICARD/MarkDuplicates.jar I=$SAMPLE_NAME.withRG.GATKRecalibrated.bam O=$SAMPLE_NAME.withRG.GATKRecalibrated.flagged.bam METRICS_FILE=$SAMPLE_NAME.Duplicates_metrics.txt VALIDATION_STRINGENCY=SILENT TMP_DIR=tmp/ 

$SAMTOOLS index $SAMPLE_NAME.withRG.GATKRecalibrated.flagged.bam
#rm -f $SAMPLE_NAME.withRG.GATKRecalibrated.bam
#rm -f $SAMPLE_NAME.withRG.GATKRecalibrated.bam.bai
#rm -f $SAMPLE_NAME.Duplicates_metrics.txt

echo 'Pre-processing complete'

echo python $PRADA_FUSION -bam $SAMPLE_NAME.withRG.GATKRecalibrated.flagged.bam -conf $PRADA_CONF_F -tag $SAMPLE_NAME -junL 80 -outdir $OUT_DIR/fusions
/usr/local/python-2.7.2/bin/python2.7 $PRADA_FUSION -bam $SAMPLE_NAME.withRG.GATKRecalibrated.flagged.bam -conf $PRADA_CONF_F -tag $SAMPLE_NAME -junL 80 -outdir $OUT_DIR/fusions
