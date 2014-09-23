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

def emit_fasta(seqs):
  for name,seq in seqs.items():
    print(">" + name)
    print(seq)

def f2p(seqs):
  num = len(seqs)
  length = len(list(seqs.values())[0])
  p_seqs = []
  for i,(name,seq) in enumerate(seqs.items()):
    p_seq = ""
    p_seq += name[:9]
    p_seq += " "
    p_seq += seq
    p_seqs.append(p_seq)
  return p_seqs, num, length

def phylnames(seqs):
  sorted_keys = sorted(seqs.keys())
  p_chain = {}
  for key in sorted_keys:
    if key[1:9] not in p_chain:
      p_chain[key[1:9]] = []
    p_chain[key[1:9]].append(key)

  p_seqs = {}
  alpha = list(string.ascii_lowercase) + list(string.ascii_uppercase)
  for name, key_list in p_chain.items():
    for i,key in enumerate(key_list):
      p_seqs[name + alpha[i]] = seqs[key]
  return p_seqs

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('fasta_file', help="the fasta file to convert")
  parser.add_argument('-f', action='store_true', help="just rename the sequences to phylip")
  args = parser.parse_args()
  seqs = parse_fasta(args.fasta_file)
  seqs = phylnames(seqs)
  if not args.f:
    p_seqs, num, length = f2p(seqs)
    print(" " + str(num) + " " + str(length))
    for line in p_seqs:
      print(line)
  else:
    emit_fasta(seqs)
