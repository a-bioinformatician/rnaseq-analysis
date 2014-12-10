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

#REFERENCE FILES
TRX_G_FA=/ifs/rcgroups/ccgd/rpa4/software/pyPRADA_1.2/PRADA-ref-hg19/Ensembl64.transcriptome.plus.genome.fasta
TRX_G_MAP=/ifs/rcgroups/ccgd/rpa4/software/pyPRADA_1.2/PRADA-ref-hg19/Ensembl64.transcriptome.plus.genome.map
HG19_FA=/ifs/rcgroups/ccgd/rpa4/software/pyPRADA_1.2/PRADA-ref-hg19/Homo_sapiens_assembly19.fasta
DBSNP_VCF=/ifs/rcgroups/ccgd/rpa4/software/pyPRADA_1.2/PRADA-ref-hg19/dbsnp_135.b37.vcf
HG19_GTF=/ifs/rcgroups/ccgd/rpa4/software/pyPRADA_1.2/PRADA-ref-hg19/Homo_sapiens.GRCh37.64.gtf

# SAMPLE INFORMATION
SAMPLE_NAME=$1
SAMPLE_BAM_PATH=$2 # PATH TO SAMPLE BAM FILE.
FASTQ_OUTDIR=$3
OUT_DIR=$4

#PIPELINE PASS IN FQ1 and FQ2

# GENERATE FASTQ DATA FROM ALIGNED SAMPLE BAM FILE
if [ -e $FASTQ_OUTDIR/'read1.fastq' ] && [ -e $FASTQ_OUTDIR/'read2.fastq' ] ; then
  echo 'Fastq files previously generated'
else 
  echo 'Generating fastq files from BAM file.'
  sh $BAM2FQ $SAMPLE_BAM_PATH $FASTQ_OUTDIR
fi

mkdir -p $OUT_DIR/temp
cd $OUT_DIR
for READ in 1 2 
do
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # REALIGN THE READS TO THE TRANSCRIPTOME AND GENOME.
  echo bwa aln -t 8 $TRX_G_FA $FASTQ_OUTDIR/'read'$READ'.fastq' '->' $OUT_DIR/temp/'read'$READ'.sai'
#  $BWA aln -t 8 $TRX_G_FA $FASTQ_OUTDIR/'read'$READ'.fastq' > $OUT_DIR/temp/'read'$READ'.sai'

  echo bwa samse -s -n 100 $TRX_G_FA $OUT_DIR/temp/'read'$READ'.sai' $FASTQ_OUTDIR/'read'$READ'.fastq' '->' $OUT_DIR/temp/'read'$READ'.sam'
#  $BWA samse -s -n 100 $TRX_G_FA $OUT_DIR/temp/'read'$READ'.sai' $FASTQ_OUTDIR/'read'$READ'.fastq' > $OUT_DIR/temp/'read'$READ'.sam'
#  rm $OUT_DIR/temp/'read'$READ'.sai'

  echo samtools view -bS -o 'read'$READ'.bam' temp/'read'$READ'.sam'
#  samtools view -bS -o 'read'$READ'.bam' temp/'read'$READ'.sam'
#  rm $OUT_DIR/temp/'read'$READ'.sam'

  echo samtools sort -n -m 1000000000 'read'$READ'.bam' 'read'$READ'.sorted'
#  samtools sort -n -m 1000000000 'read'$READ'.bam' 'read'$READ'.sorted'
#  rm 'read'$READ'.bam'
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # REMAP WITH GATK
#  java -Djava.io.tmpdir=tmp/ -cp $GATK/RemapAlignments.jar -Xmx8g org.broadinstitute.cga.tools.gatk.rna.RemapAlignments M=$TRX_G_MAP IN='read'$READ'.sorted.bam' OUT='read'$READ'.remapped.bam' R=$HG19_FA REDUCE=TRUE
#  rm -f 'read'$READ'.sorted.bam'

#  samtools sort -n -m 1000000000 'read'$READ'.remapped.bam' 'read'$READ'.remapped.sorted'
#  rm -f 'read'$READ'.remapped.bam'
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REJOIN PAIRED READS INTO SINGLE BAM FILE
#java -Djava.io.tmpdir=tmp/ -Xmx8g -jar $GATK/PairMaker.jar IN1='read1.remapped.sorted.bam' IN2='read2.remapped.sorted.bam' OUTPUT=$SAMPLE_NAME.paired.bam TMP_DIR=tmp/
#rm -f 'read1.remapped.sorted.bam'
#rm -f 'read2.remapped.sorted.bam'

#samtools sort -m 1000000000 $SAMPLE_NAME.paired.bam $SAMPLE_NAME.paired.sorted
#rm -f $SAMPLE_NAME.paired.bam
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#java -Xmx8g -jar $PICARD/AddOrReplaceReadGroups.jar I=$SAMPLE_NAME.paired.sorted.bam O=$SAMPLE_NAME.withRG.paired.sorted.bam RGLB=$SAMPLE_NAME RGPL=Illumina RGPU=$SAMPLE_NAME RGSM=$SAMPLE_NAME
#rm -f $SAMPLE_NAME.paired.sorted.bam

#samtools index $SAMPLE_NAME.withRG.paired.sorted.bam
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar -l INFO -R $HG19_FA --default_platform illumina --knownSites $DBSNP_VCF -I $SAMPLE_NAME.withRG.paired.sorted.bam --downsample_to_coverage 10000 -T CountCovariates -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -nt 8 -recalFile $SAMPLE_NAME.orig.csv

#java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar -l INFO -R $HG19_FA --default_platform illumina -I $SAMPLE_NAME.withRG.paired.sorted.bam -T TableRecalibration --out $SAMPLE_NAME.withRG.GATKRecalibrated.bam -recalFile $SAMPLE_NAME.orig.csv
#rm -f $SAMPLE_NAME.withRG.paired.sorted.bam
#rm -f $SAMPLE_NAME.withRG.paired.sorted.bam.bai
#rm -f $SAMPLE_NAME.orig.csv

#java -Xmx8g -jar $PICARD/MarkDuplicates.jar I=$SAMPLE_NAME.withRG.GATKRecalibrated.bam O=$SAMPLE_NAME.withRG.GATKRecalibrated.flagged.bam METRICS_FILE=$SAMPLE_NAME.Duplicates_metrics.txt VALIDATION_STRINGENCY=SILENT TMP_DIR=tmp/ 

#samtools index $SAMPLE_NAME.withRG.GATKRecalibrated.flagged.bam
#rm -f $SAMPLE_NAME.withRG.GATKRecalibrated.bam
#rm -f $SAMPLE_NAME.withRG.GATKRecalibrated.bam.bai
#rm -f $SAMPLE_NAME.Duplicates_metrics.txt

gtf=/ifs/rcgroups/ccgd/rpa4/analysis/rnaseq/SRS0001_PDL1-Therapy-LungCA/Run287_Lane1_D381plus_BC1/prada/test.gtf
# gtf=/ifs/rcgroups/ccgd/rpa4//ref2/human/gencode/current/annotation/raw/gencode.v17.annotation.gtf
# gtf=/ifs/rcgroups/ccgd/rpa4/analysis/rnaseq/SRS0001_PDL1-Therapy-LungCA/Run287_Lane1_D381plus_BC1/prada/prada_sorted.gtf
# fa=/ifs/rcgroups/ccgd/rpa4//ref2/human/gencode/current/genome/primary-assembly/gencode_GRCh37-p11_dna.primary-assembly.fa
# fa=/ifs/rcgroups/ccgd/rpa4//ref2/human/ensembl/grch37-75/genome/primary-assembly/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
# fa=$HG19_FA
fa=/ifs/rcgroups/ccgd/rpa4/software/pyPRADA_1.2/PRADA-ref-hg19/ensembl_test.ordered.fasta

# echo java -jar $PICARD/ReorderSam.jar INPUT=$TEMP_DIR/alignments.rg.bam OUTPUT=$TEMP_DIR/alignments.rg.reorder.bam REFERENCE=$REF_FASTA_BASE.fa
# java -jar $PICARD/ReorderSam.jar INPUT=$SAMPLE_NAME.withRG.GATKRecalibrated.flagged.bam OUTPUT=$SAMPLE_NAME.withRG.GATKRecalibrated.flagged.reordered.bam REFERENCE=$fa

# echo java -jar $PICARD/BuildBamIndex.jar INPUT=$SAMPLE_NAME.withRG.GATKRecalibrated.flagged.reordered.bam
java -jar $PICARD/BuildBamIndex.jar INPUT=$SAMPLE_NAME.withRG.GATKRecalibrated.flagged.reordered.bam 

echo java -Xmx8g -jar $QC -ttype 2 -t $HG19_GTF -r $HG19_FA -s "$SAMPLE_NAME|$SAMPLE_NAME.withRG.GATKRecalibrated.flagged.reordered.bam|$SAMPLE_NAME" -o $SAMPLE_NAME -n 1000
java -Xmx8g -jar $QC -t $gtf -r $fa -s "$SAMPLE_NAME|$SAMPLE_NAME.withRG.GATKRecalibrated.flagged.reordered.bam|$SAMPLE_NAME" -o $SAMPLE_NAME -n 1000

echo 'Pre-processing complete'

echo python $PRADA_FUSION -bam test_sample.withRG.GATKRecalibrated.flagged.bam -conf /ifs/rcgroups/ccgd/software/pyPRADA_1.2/conf.txt -tag test_sample -junL 80 -outdir $PWD
/usr/local/python-2.7.2/bin/python2.7 $PRADA_FUSION -bam test_sample.withRG.GATKRecalibrated.flagged.bam -conf /ifs/rcgroups/ccgd/software/pyPRADA_1.2/conf.txt -tag test_sample -junL 80 -outdir $PWD
