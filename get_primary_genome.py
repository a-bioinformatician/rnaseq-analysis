#! /usr/bin/local/python

'''
Author: Ryan Abo
Date: 12/15/2014

This script is designed to take in a single fasta file and only keep the primary chromosomes
plus the unplaced/unlocalized scaffolds. 

It is mainly designed for the Gencode genome builds that include all the top-level chromosomes
and unplaced/unlocalized as well as patches and haplotypes and alternative chroms.
(ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/)

The patches and alternate assemblies/haplotypes are NOT kept.

INPUT: 
    1. fasta file with multiple seq records.
    2. text file with a list of chromosomes to keep.
OUTPUT: 
    1. fasta file with primary seq records.
'''

import sys, os
from Bio import SeqIO
from optparse import OptionParser

usage = '%prog [options] <fasta file name> <text file with chrom record ids to keep> <output fasta file name>'
desc = """Script to parse fasta file and keep primary assembly chroms."""
parser = OptionParser(usage=usage,description=desc)

def write_fa(rec) :
  print rec.id
  fname = rec.id
  if fname.find('chr') == -1 : fname = 'chr' + fname
  fout = open(fname+".fa", 'w')
  SeqIO.write(rec,fout,"fasta")
  fout.close()

if __name__ == '__main__' :
    opts, args = parser.parse_args(sys.argv[1:])
    fasta_fn = args[0]
    chrom_lst_fn = args[1]
    out_fasta_fn = args[2]

    chrom_keep_lst = []
    for line in open(chrom_lst_fn, 'rU').readlines() :
        line = line.strip()
        chrom_keep_lst.append(line)

    fasta_f = open(fasta_fn, 'rU')
    new_fasta_f = open(out_fasta_fn, 'w')
    keep_recs = []
    for rec in SeqIO.parse(fasta_f,"fasta") :
        rec.id = rec.id.split()[0]
        if rec.id in chrom_keep_lst :
            keep_recs.append(rec)
    SeqIO.write(keep_recs, new_fasta_f, 'fasta')
    fasta_f.close()
    new_fasta_f.close()
