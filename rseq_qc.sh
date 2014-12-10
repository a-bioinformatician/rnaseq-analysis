#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M ryanp_abo@dfci.harvard.edu
#$ -pe pvm 8
###################

# ------------------------------
# Script to run RSeQC on alignments 
# ------------------------------

# VARIABLES
SAMPLE_NAME=$1
BAM=$2
OUT_DIR=$3

# PATHS

REF_PATH='/ifs/rcgroups/ccgd/rpa4/ref2/human/gencode'
REF_NAME=`cat $REF_PATH/current_version`
REF_BASE=$REF_PATH/current/genome/primary-assembly/$REF_NAME
RNA_BED='/ifs/rcgroups/ccgd/rpa4/ref2/human/ucsc/hg19_rRNA.bed'
GTF_BED=$REF_PATH/current/annotation/primary-assembly/$REF_NAME.bed

mkdir -p $OUT_DIR
cd $OUT_DIR

echo bam_stat.py -i $BAM
bam_stat.py -i $BAM >> bam_stats.log 2>&1

echo clipping_profile.py -i $BAM -o $SAMPLE_NAME
clipping_profile.py -i $BAM -o $SAMPLE_NAME

echo geneBody_coverage.py -r $GTF_BED -i $BAM -o $SAMPLE_NAME
geneBody_coverage.py -r $GTF_BED -i $BAM -o $SAMPLE_NAME

echo infer_experiment.py -r $GTF_BED -i $BAM
infer_experiment.py -r $GTF_BED -i $BAM >> infer_experiment.log 2>&1

echo inner_distance.py -i $BAM -o $SAMPLE_NAME -r $GTF_BED 
inner_distance.py -i $BAM -o $SAMPLE_NAME -r $GTF_BED 

echo junction_annotation.py -i $BAM -r $GTF_BED -o $SAMPLE_NAME
junction_annotation.py -i $BAM -r $GTF_BED -o $SAMPLE_NAME >> junction_annotation.log 2>&1

echo junction_saturation.py -i $BAM -o $SAMPLE_NAME -r $GTF_BED  
junction_saturation.py -i $BAM -o $SAMPLE_NAME -r $GTF_BED  

echo read_distribution.py -i $BAM -r $GTF_BED
read_distribution.py -i $BAM -r $GTF_BED >> read_distribution.log 2>&1

echo read_duplication.py -i $BAM -o $SAMPLE_NAME
read_duplication.py -i $BAM -o $SAMPLE_NAME

echo read_GC.py -i $BAM -o $SAMPLE_NAME
read_GC.py -i $BAM -o $SAMPLE_NAME

echo read_NVC.py -i $BAM -o $SAMPLE_NAME
read_NVC.py -i $BAM -o $SAMPLE_NAME

echo read_quality.py -i $BAM -o $SAMPLE_NAME
read_quality.py -i $BAM -o $SAMPLE_NAME

echo read_hexamer.py -i $BAM -r $GTF_BED
read_hexamer.py -i $BAM -r $GTF_BED

echo split_bam.py -i $BAM -r $RNA_BED -o $SAMPLE_NAME
split_bam.py -i $BAM -r $RNA_BED -o $SAMPLE_NAME >> split_bam.log 2>&1
