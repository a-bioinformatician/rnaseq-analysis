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
# SAMPLE_NAME=$1
# SAMPLE_BAM_PATH=$2 # PATH TO SAMPLE BAM FILE.
# FASTQ_OUTDIR=$3
# OUT_DIR=$4

python $PRADA_FUSION -bam test_sample.withRG.GATKRecalibrated.flagged.bam -conf /ifs/rcgroups/ccgd/software/pyPRADA_1.2/conf.txt -tag test_sample -junL 80 -outdir $PWD
