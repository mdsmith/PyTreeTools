#! /usr/bin/env python3

import argparse

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

def f2p(seqs):
  num = len(seqs)
  length = len(list(seqs.values())[0])
  p_seqs = []
  for i,(name,seq) in enumerate(seqs.items()):
    p_seq = ""
    p_seq += name[1:10]
    p_seq += " "
    p_seq += seq
    p_seqs.append(p_seq)
  return p_seqs, num, length

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('fasta_file', help="the fasta file to convert")
  args = parser.parse_args()
  seqs = parse_fasta(args.fasta_file)
  p_seqs, num, length = f2p(seqs)
  print(" " + str(num) + " " + str(length))
  for line in p_seqs:
    print(line)
