#!/usr/bin/env python

import sys
import os

if __name__ == "__main__" :
  if len(sys.argv) != 3:
    print "Parses STAR alignment log.\n"
    print "usage:", sys.argv[0], "<IN FILE> <OUT FILE>"
    sys.exit()
 
  infile = sys.argv[1]
  outfile = sys.argv[2]

  f = open(infile, 'r')
  o = open(outfile, 'w')

  hit = 0
  cat = ''
  for line in f.readlines() :
    line = line.strip()
    linesplit = line.split()
    if line.find('Mapping speed') > -1 : 
      hit = 1
      continue

    if hit == 1 : hit += 1
    elif hit > 1 :
      if line.find('UNIQUE') > -1 : 
        cat = 'Unique reads'
      elif line.find('MULTI-MAPPING') > -1 : 
        cat = 'Multi-mapped reads'
      elif line.find('UNMAPPED READS') > -1 :
        cat = 'Unmapped reads'
      else :
        lsplit = line.split('|')
        out = [cat, lsplit[0], lsplit[1].lstrip('\t')]
        o.write('\t'.join(out) + '\n')
      hit += 1     

  o.close()
  f.close()
