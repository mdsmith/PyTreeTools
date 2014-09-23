#! /usr/bin/env python3

import argparse
import re

def parse_fasta(fasta_file_name):
  fh = open(fasta_file_name, 'r')
  lines = fh.readlines()
  return lines

def clean_seqs(seqs):
  for i,line in enumerate(seqs):
    if line[0] == '>':
      line = re.sub('[^a-zA-Z0-9\n\>]', '_', line)
    line = line.strip('\n')
    seqs[i] = line
  return seqs

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('fasta_file', help="the fasta file to clean")
  args = parser.parse_args()
  seqs = parse_fasta(args.fasta_file)
  clean = clean_seqs(seqs)
  for line in clean:
    print(line)
  #print(clean)
