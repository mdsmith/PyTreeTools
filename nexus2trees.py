#! /usr/bin/env python3

import argparse
from parsers import parse_nexus

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('nexus_file', help="the nexus file to convert")
  args = parser.parse_args()
  seqs, trees = parse_nexus(args.nexus_file)
  if not len(seqs.items()):
    print("No sequences found")
  for name, newick in trees.items():
    print(newick.strip().replace(';', ''))

