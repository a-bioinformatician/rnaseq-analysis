#!/usr/bin/env python

import sys
import os
import glob

if __name__ == "__main__" :
  if len(sys.argv) != 5:
    print "Parses sample summary files and merges them together.\n"
    print "usage:", sys.argv[0], "<PATH> <SAMPLE NAMES> <OUT DIR> <OUT NAME>"
    sys.exit()
 
  analysis_dir = sys.argv[1]
  sample_names = sys.argv[2].split(',')
  out_dir = sys.argv[3]
  out_name = sys.argv[4]

  # RAW FLAGSTAT SUMMARIES
  flagstat_out = open(os.path.join(out_dir, out_name + '.flagstat_summary.txt'), 'w')
  header = ['Raw sequence summary']
  values = []
  iter = 0
  for sample in sample_names : 
    fnames = glob.glob(os.path.join(analysis_dir, sample, 'summary', '*.bam.flagstat_summary.txt'))
    for fname in fnames :
      sample_fname = os.path.basename(fname).split('.')[0] 
      f = open(fname, 'rU')
      line_iter = 0
      for line in f.readlines() :
        line = line.strip()
        lsplit = line.split('\t')
        if iter == 0 : 
          values.append(lsplit)
        else :
          values[line_iter].append(lsplit[1])
        line_iter += 1
      header.append(sample_fname)
      iter += 1  
  flagstat_out.write('\t'.join(header) + '\n')
  for value in values :
    flagstat_out.write('\t'.join(value) + '\n')
  flagstat_out.close()

  # EXPERIMENT SUMMARIES
  experiment_out = open(os.path.join(out_dir, out_name + '.experiment.txt'), 'w')
  header = ['RSEQC experiment type']
  values = []
  iter = 0
  for sample in sample_names :
    fname = os.path.join(analysis_dir, sample, 'summary', sample + '.experiment.txt')
    print fname
    sample_fname = os.path.basename(fname).split('.')[0]
    f = open(fname, 'rU')
    line_iter = 0
    for line in f.readlines() :
      line = line.strip()
      print line
      lsplit = line.split('\t')
      if iter == 0 :
        values.append(lsplit)
      else :
        values[line_iter].append(lsplit[1])
      line_iter += 1
    header.append(sample_fname)
    iter += 1
  experiment_out.write('\t'.join(header) + '\n')
  for value in values :
    experiment_out.write('\t'.join(value) + '\n')
  experiment_out.close()

  # TRX ALIGN SUMMARIES
  fout = open(os.path.join(out_dir, out_name + '.trx_align_summary.txt'), 'w')
  header = ['Transcriptome alignment']
  values = []
  iter = 0
  for sample in sample_names :
    fname = os.path.join(analysis_dir, sample, 'summary', sample + '.trx_align_summary.txt')
    sample_fname = os.path.basename(fname).split('.')[0]
    f = open(fname, 'rU')
    line_iter = 0
    for line in f.readlines() :
      line = line.strip()
      lsplit = line.split('\t')
      if iter == 0 :
        values.append(lsplit)
      else :
        values[line_iter].append(lsplit[1])
      line_iter += 1
    header.append(sample_fname)
    iter += 1
  fout.write('\t'.join(header) + '\n')
  for value in values :
    fout.write('\t'.join(value) + '\n')
  fout.close()

  # READ DISTRIBUTION SUMMARIES
  rd_out = open(os.path.join(out_dir, out_name + '.read_distribution.txt'), 'w')
  header = ['RSEQC read distribution']
  values = []
  iter = 0
  for sample in sample_names :
    fname = os.path.join(analysis_dir, sample, 'summary', sample + '.read_distribution.txt')
    sample_fname = os.path.basename(fname).split('.')[0]
    f = open(fname, 'rU')
    line_iter = 0
    for line in f.readlines() :
      line = line.strip()
      lsplit = line.split('\t')
      if iter == 0 :
        values.append(lsplit[0:2])
      else :
        values[line_iter].append(lsplit[1])
      line_iter += 1
    header.append(sample_fname)
    iter += 1
  rd_out.write('\t'.join(header) + '\n')
  for value in values :
    rd_out.write('\t'.join(value) + '\n')
  rd_out.close()

  # RIBOSOMAL RNA SUMMARIES 
  fout = open(os.path.join(out_dir, out_name + '.rRNA_alignment.txt'), 'w')
  header = ['rRNA alignment %']
  values = []
  iter = 0
  for sample in sample_names :
    fname = os.path.join(analysis_dir, sample, 'summary', sample + '.rRNA_align.flagstat_summary.txt')
    sample_fname = os.path.basename(fname).split('.')[0]
    f = open(fname, 'rU')
    line_iter = 0
    for line in f.readlines() :
      line = line.strip()
      lsplit = line.split('\t')
      if iter == 0 :
        values.append(lsplit)
      else :
        values[line_iter].append(lsplit[1])
      line_iter += 1
    header.append(sample_fname)
    iter += 1
  fout.write('\t'.join(header) + '\n')
  for value in values :
    fout.write('\t'.join(value) + '\n')
  fout.close()

  # STAR ALiGNMENT SUMMARIES 
  fout = open(os.path.join(out_dir, out_name + '.star_alignment.txt'), 'w')
  header = ['STAR alignment','']
  values = []
  iter = 0
  for sample in sample_names :
    fname = os.path.join(analysis_dir, sample, 'summary', sample + '.star_log_summary.txt')
    sample_fname = os.path.basename(fname).split('.')[0]
    f = open(fname, 'rU')
    line_iter = 0
    for line in f.readlines() :
      line = line.strip()
      lsplit = line.split('\t')
      if len(lsplit) == 2 : lsplit.insert(0,'')
      if iter == 0 :
        if lsplit[2].find('%') > -1 and lsplit[1].find('%') == -1 :
          lsplit[1] += '(%)'
        lsplit[2] = lsplit[2].strip('%')
        values.append(lsplit)
      else :
        lsplit[-1] = lsplit[-1].strip('%')
        values[line_iter].append(lsplit[-1])
      line_iter += 1
    header.append(sample_fname)
    iter += 1
  fout.write('\t'.join(header) + '\n')
  for value in values :
    fout.write('\t'.join(value) + '\n')
  fout.close()
