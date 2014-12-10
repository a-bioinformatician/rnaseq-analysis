#!/usr/bin/env python

import sys
import os

if __name__ == "__main__" :
  if len(sys.argv) != 4:
    print "Parses samtools flagstat output file.\n"
    print "usage:", sys.argv[0], "<IN FILE> <OUT FILE> <FUNCTION>"
    sys.exit()
 
  infile = sys.argv[1]
  outfile = sys.argv[2]
  type = sys.argv[3]

  f = open(infile, 'r')
  o = open(outfile, 'w')
  for line in f.readlines() :
    line = line.strip()
    if type == 'rRNA' and line.find('mapped (') > -1 :
      rrna_perc = line.split('(')[1].split(':')[0].strip('%')
      o.write('rRNA alignment\t' + rrna_perc)
    elif type == 'qc' and line.find('QC-passed') > -1 : 
      qc_passed = line.split(' + ')[0].strip()
      qc_failed = line.split(' + ')[1].split()[0]
      total_reads = int(qc_passed) + int(qc_failed)
      o.write('total_reads\t' + str(total_reads) + '\nqc_passed\t' + qc_passed + '\nqc_failed\t' + qc_failed + '\n') 
  o.close()
  f.close()
