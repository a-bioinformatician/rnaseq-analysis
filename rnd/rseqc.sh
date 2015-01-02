#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M ryanp_abo@dfci.harvard.edu
#$ -pe pvm 5
###################

# ------------------------------
# Script to run RSeQC on alignments 
# ------------------------------

rseqc=/ifs/rcgroups/ccgd/rpa4/software/RSeQC/usr/local/python-2.7.2/bin #/ifs/rcgroups/ccgd/rpa4/software/
python=/usr/local/python-2.7.2/bin/python2.7

OUT_BASENAME=$1
BAM=$2
OUT_DIR=$3
GTF_BED=$4 # $REF_PATH/current/annotation/primary-assembly/$REF_NAME.bed
RNA_BED=/ifs/rcgroups/ccgd/rpa4/ref2/human/ucsc/hg19_rRNA.bed

mkdir -p $OUT_DIR
cd $OUT_DIR

echo bam_stat.py -i $BAM
$python $rseqc/bam_stat.py -i $BAM >> bam_stats.log 2>&1

echo clipping_profile.py -i $BAM -o $OUT_BASENAME
$python $rseqc/clipping_profile.py -i $BAM -o $OUT_BASENAME

echo geneBody_coverage.py -r $GTF_BED -i $BAM -o $OUT_BASENAME
$python $rseqc/geneBody_coverage.py -r $GTF_BED -i $BAM -o $OUT_BASENAME

echo infer_experiment.py -r $GTF_BED -i $BAM
$python $rseqc/infer_experiment.py -r $GTF_BED -i $BAM >> infer_experiment.log 2>&1

echo inner_distance.py -i $BAM -o $OUT_BASENAME -r $GTF_BED 
$python $rseqc/inner_distance.py -i $BAM -o $OUT_BASENAME -r $GTF_BED 

echo junction_annotation.py -i $BAM -r $GTF_BED -o $OUT_BASENAME
$python $rseqc/junction_annotation.py -i $BAM -r $GTF_BED -o $OUT_BASENAME >> junction_annotation.log 2>&1

echo junction_saturation.py -i $BAM -o $OUT_BASENAME -r $GTF_BED  
$python $rseqc/junction_saturation.py -i $BAM -o $OUT_BASENAME -r $GTF_BED  

echo read_distribution.py -i $BAM -r $GTF_BED
$python $rseqc/read_distribution.py -i $BAM -r $GTF_BED >> read_distribution.log 2>&1

echo read_duplication.py -i $BAM -o $OUT_BASENAME
$python $rseqc/read_duplication.py -i $BAM -o $OUT_BASENAME

echo read_GC.py -i $BAM -o $OUT_BASENAME
$python $rseqc/read_GC.py -i $BAM -o $OUT_BASENAME

echo read_NVC.py -i $BAM -o $OUT_BASENAME
$python $rseqc/read_NVC.py -i $BAM -o $OUT_BASENAME

echo read_quality.py -i $BAM -o $OUT_BASENAME
$python $rseqc/read_quality.py -i $BAM -o $OUT_BASENAME

echo read_hexamer.py -i $BAM -r $GTF_BED
$python $rseqc/read_hexamer.py -i $BAM -r $GTF_BED

echo split_bam.py -i $BAM -r $RNA_BED -o $OUT_BASENAME
$python $rseqc/split_bam.py -i $BAM -r $RNA_BED -o $OUT_BASENAME >> split_bam.log 2>&1
