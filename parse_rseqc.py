#!/usr/bin/env python

import sys
import os

def parse_experiment(f, o) :
  hit = False
  for line in f.readlines() :
    line = line.strip()
    if not hit :
      if line.find('Fraction') > -1 :
        hit = True
    
    if hit : 
      o.write('\t'.join(line.split(':'))+'\n')

def parse_features(f, o) :
  hit = 0
  total = 1
  for line in f.readlines() :
    line = line.strip()
    if hit > 0 :
      if hit == 1 :
        total = int(line.split()[-1])
        o.write('\t'.join(['Total assigned tags', str(total),'1.0']) + '\n')
        hit += 1
      elif hit == 2 or hit == 3 : 
        hit += 1
      elif hit == 14 :
        print 'Done'
      else : 
        lsplit = line.split()
        out = [lsplit[0], lsplit[2], str(round((float(lsplit[2])/float(total)) * 100, 2))]
        o.write('\t'.join(out) + '\n')
        hit += 1
    else :
      if line.find('Total Tags') > -1 :
        hit = 1

if __name__ == "__main__" :
  if len(sys.argv) != 3:
    print "Parses RSEQC output file.\n"
    print "usage:", sys.argv[0], "<IN FILE> <OUT FILE> <FUNCTION>"
    sys.exit()
 
  infile = sys.argv[1]
  outfile = sys.argv[2]

  f = open(infile, 'r')
  o = open(outfile, 'w')
  if infile.find('infer_experiment') > -1 :
    parse_experiment(f, o)
  elif infile.find('read_distribution') > -1 :
    parse_features(f, o)

  o.close()
  f.close()
