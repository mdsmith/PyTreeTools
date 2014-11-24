#! /usr/bin/env python3
#
# Merge (concatenate) sequences under a single sequence name
# Particularly useful for merging two genes.
# Fasta names must match.

import argparse

def join_fastas(file1, file2):
  fh1 = open(file1, 'r')
  fh2 = open(file2, 'r')
  fh1_lines = fh1.readlines()
  fh2_lines = fh2.readlines()
  out = ""
  name = fh1_lines[0]
  seq = ""
  for i,line in enumerate(fh1_lines[1:]):
    if i == len(fh1_lines) - 2:
      seq += line
    if line[0] == '>' or i == len(fh1_lines) - 2:
      seq2 = ""
      on = False
      for line2 in fh2_lines:
        if not on:
          if line2 == name:
            on = True
        else:
          if line2[0] == '>':
            break
          else:
            seq2 += line2
      out += name
      out += seq
      out += seq2
      name = line
      seq = ""
    else:
      seq += line
  return out


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('file1', help="the first file")
  parser.add_argument('file2', help="the second file")
  args = parser.parse_args()

  print(join_fastas(args.file1, args.file2))
