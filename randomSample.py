#! /usr/bin/env python3

import argparse
import sys
import random

def parse_fasta(fasta_file_name):
  fh = open(fasta_file_name, 'r')
  lines = fh.readlines()
  seqs = {}
  seq = ""
  name = ""
  for i,line in enumerate(lines):
    if line[0] == '>' or i == len(lines) - 1:
      if i == len(lines) - 1:
        seq += line
      if i != 0:
        seqs[name.strip('\n')] = seq.strip('\n')
      name = ""
      seq = ""
      name = line
    else:
      seq += line
  return seqs

def sample_fasta(seqs, number):
  to_pop = len(seqs) - int(number)
  for _ in range(to_pop):
    seqs.pop(random.choice(list(seqs.keys())))
  return seqs

def print_dict_to_fasta(seqs):
  for key, value in seqs.items():
    print(key)
    print(value)
    print('')

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('fasta_file', help="the fasta file to clean")
  parser.add_argument('number', help="the number of sequences to sample \
      (without replacement)")
  parser.add_argument('-d', action='store_true', help="write dupes to stderr")
  args = parser.parse_args()
  seqs = parse_fasta(args.fasta_file)
  sampled = sample_fasta(seqs, args.number)
  print_dict_to_fasta(sampled)
