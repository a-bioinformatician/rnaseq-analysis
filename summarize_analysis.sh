#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -pe pvm 4
###################

# ------------------------------
# ------------------------------


DIR=$1
OUT_DIR=$2
SAMPLE_NAME=$3

# SOFTWARE
samtools='/usr/local/samtools-0.1.18-default/samtools'
parse_flagstat='/ifs/rcgroups/ccgd/rpa4/scripts/github_repos/rnaseq-analysis/parse_flagstat.py'
parse_rseqc='/ifs/rcgroups/ccgd/rpa4/scripts/github_repos/rnaseq-analysis/parse_rseqc.py'
parse_star_log='/ifs/rcgroups/ccgd/rpa4/scripts/github_repos/rnaseq-analysis/parse_star_log.py'
parse_rsem_log='/ifs/rcgroups/ccgd/rpa4/scripts/github_repos/rnaseq-analysis/parse_rsem_log.py'

# PARAMS
mkdir -p $OUT_DIR

BAMS=`ls $DIR/data | grep .bam`
echo 'BAMs' $BAMS
for BAM in $BAMS;
do
  # Unaligned reads stats
  echo samtools flagstat $DIR/data/$BAM '>' $OUT_DIR/$BAM.flagstat.txt
  samtools flagstat $DIR/data/$BAM > $OUT_DIR/$BAM.flagstat.txt
  python $parse_flagstat $OUT_DIR/$BAM.flagstat.txt $OUT_DIR/$BAM.flagstat_summary.txt qc
done

# Genome alignment stats
echo python $parse_star_log $DIR/star/Log.final.out $OUT_DIR/star_log_summary.txt
python $parse_star_log $DIR/star/Log.final.out $OUT_DIR/$SAMPLE_NAME.star_log_summary.txt

echo samtools flagstat $DIR/rRNA_align/rRNA_align.sorted.bam '>' $OUT_DIR/$SAMPLE_NAME.rRNA_align.flagstat.txt
samtools flagstat $DIR/rRNA_align/hgRibo_align.sorted.bam > $OUT_DIR/$SAMPLE_NAME.rRNA_align.flagstat.txt
echo python $parse_flagstat $OUT_DIR/$SAMPLE_NAME.flagstat.txt $OUT_DIR/$SAMPLE_NAME.rRNA_align.flagstat_summary.txt rRNA
python $parse_flagstat $OUT_DIR/$SAMPLE_NAME.rRNA_align.flagstat.txt $OUT_DIR/$SAMPLE_NAME.rRNA_align.flagstat_summary.txt rRNA

echo python $parse_rseqc $DIR/qc/rseqc/infer_experiment.log $OUT_DIR/$SAMPLE_NAME.experiment.txt
python $parse_rseqc $DIR/qc/rseqc/infer_experiment.log $OUT_DIR/$SAMPLE_NAME.experiment.txt
echo python $parse_rseqc $DIR/qc/rseqc/read_distribution.log $OUT_DIR/$SAMPLE_NAME.read_distribution.txt
python $parse_rseqc $DIR/qc/rseqc/read_distribution.log $OUT_DIR/$SAMPLE_NAME.read_distribution.txt

# Transcriptome alignment stats
echo 'Extracting counts from transcriptome alignment'
python $parse_rsem_log $DIR/logs/rsem_expression.sh.e* $OUT_DIR/$SAMPLE_NAME.trx_align_summary.txt


