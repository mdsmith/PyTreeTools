#! /usr/bin/env python3

import argparse
import sys

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

def rem_dupes(seqs, to_err):
  keys_to_pop = []
  for i, (key1, value1) in enumerate(seqs.items()):
    for key2, value2 in seqs.items():
      if value1 == value2 and key1 is not key2:
        keys_to_pop.append(key2)
        if to_err:
          print(key1, file=sys.stderr)
          print(value1, file=sys.stderr)
          print(key2, file=sys.stderr)
          print(value2, file=sys.stderr)
  for key in keys_to_pop:
    seqs.pop(key, None)
  return seqs

def print_dict_to_fasta(seqs):
  for key, value in seqs.items():
    print(key)
    print(value)
    print('')

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('fasta_file', help="the fasta file to clean")
  parser.add_argument('-d', action='store_true', help="write dupes to stderr")
  args = parser.parse_args()
  seqs = parse_fasta(args.fasta_file)
  clean = rem_dupes(seqs, args.d)
  print_dict_to_fasta(clean)
