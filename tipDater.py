#! /usr/bin/env python3

import argparse
import sys

def parse_file(fasta_file_name):
  fh = open(fasta_file_name, 'r')
  lines = fh.readlines()
  return lines

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

def parse_csv(csv_file_name, header):
  fh = open(csv_file_name, 'r')
  lines = fh.readlines()
  if header:
    lines = lines[1:]
  dates = {}
  for line in lines:
    #print(line)
    try:
      key, value = line.split(',')
      dates[key] = value.strip()
    except ValueError:
      print("broken date: " + line, file=sys.stderr)

  return dates

def date_tips(seqs, dates, omit):
  dated = {}
  for name, seq in seqs.items():
    seq_found = False
    for accn, date in dates.items():
      if accn in name:
        dated_name = name + "_" + date
        dated[dated_name] = seq
        seq_found = True
    if seq_found == False:
      print("not dated: " + name, file=sys.stderr)
      if not omit:
        dated[name+"_x"] = seq
  return dated

def print_dict_to_fasta(seqs):
  for key, value in seqs.items():
    print(key)
    print(value)
    print('')

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('infile', help="the file to date")
  parser.add_argument('datefile', help="the csv containing name, csv")
  parser.add_argument('-c',
                      dest = "header",
                      action='store_true',
                      help="ignore any csv header")
  parser.add_argument('-o',
                      dest = "omit",
                      action='store_true',
                      help="omit undated sequences")
  args = parser.parse_args()
  lines = parse_fasta(args.infile)
  dates = parse_csv(args.datefile, args.header)
  dated = date_tips(lines, dates, args.omit)
  print_dict_to_fasta(dated)
