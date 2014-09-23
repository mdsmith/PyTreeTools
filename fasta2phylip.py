#! /usr/bin/env python3

import argparse
import sys
import string

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
    p_seq += name[1:9]
    p_seq += "  "
    p_seq += seq
    p_seqs.append(p_seq)
  return p_seqs, num, length

def check_dupes(seqs):
  alpha = list(string.ascii_lowercase) + list(string.ascii_uppercase)
  for i,seq in enumerate(seqs):
    found = 0
    name1 = seq[:10]
    for i2,seq2 in enumerate(seqs[i+1:]):
      name2 = seq2[:10]
      if name1 == name2:
        found += 1
        print("dupe " + str(found) + " found: " + name1 + " " + name2, file=sys.stderr)
        print(seqs[i2 + i + 1][:12], file=sys.stderr)
        seqs[i2 + i + 1] = seqs[i2 + i + 1][:8] + alpha[found] + seqs[i2 + i
            + 1][9:]
        print(seqs[i2 + i + 1][:12], file=sys.stderr)


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('fasta_file', help="the fasta file to convert")
  args = parser.parse_args()
  seqs = parse_fasta(args.fasta_file)
  p_seqs, num, length = f2p(seqs)
  check_dupes(p_seqs)
  print(" " + str(num) + " " + str(length))
  for line in p_seqs:
    print(line)
