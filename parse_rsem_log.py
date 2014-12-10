#!/usr/bin/env python

import sys
import os

if __name__ == "__main__" :
  if len(sys.argv) != 3:
    print "Parses RSEM Bowtie alignment log file.\n"
    print "usage:", sys.argv[0], "<IN FILE> <OUT FILE>"
    sys.exit()
 
  infile = sys.argv[1]
  outfile = sys.argv[2]

  f = open(infile, 'r')
  o = open(outfile, 'w')
  for line in f.readlines() :
    line = line.strip()
    if line.find("#") > -1 : 
      line = line.split('# ')[1]
      lsplit = line.split(": ")
      if lsplit[0].find('reads processed') > -1 :
        o.write('Number of reads processed\t' + lsplit[1] + '\n')
      elif lsplit[0].find('reads with at least one') > -1 :
        o.write('Number of reads with at least one alignment to trx\t' + lsplit[1].split()[0] + '\n')
      elif lsplit[0].find('reads that failed') > -1 : 
        o.write('Number of reads with no alignment to trx\t' + lsplit[1].split()[0] + '\n')
        break  
  o.close()
  f.close()
