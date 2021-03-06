#! /usr/bin/env python3

import argparse
import re
import parsers

def parse_file(fasta_file_name):
  fh = open(fasta_file_name, 'r')
  lines = fh.readlines()
  return lines

def clean_seqs(contents, name, num_ok=False):
  if contents[0][0] == ">":
    return clean_fasta(contents)
  else:
    return clean_newick(contents, name)

def clean_fasta(seqs):
  for i,line in enumerate(seqs):
    if line[0] == '>':
      line = clean_name(line)
      #line = re.sub('[^a-zA-Z0-9\n\>]', '_', line)
    line = line.strip('\n')
    seqs[i] = line
  return seqs

def clean_newick(contents, name):
  in_tree = parsers.Tree(contents[0])
  if name:
    clean_node_name(in_tree.root, 1)
  else:
    clean_node_name(in_tree.root, -1)
  return in_tree.to_newick_string()

def clean_node_name(node, name):
  node.name = clean_name(node.name)
  if name != -1:
    if node.name == "":
      node.name = "Node" + str(name)
      name += 1
  for c in node.children:
    clean_node_name(c, name)

def is_num(s):
  try:
    float(s)
    return True
  except ValueError:
    return False

def clean_name(inname, num_ok=False):
  if is_num(inname) and not num_ok:
    return ""
  elif is_num(inname) and num_ok:
    return inname
  else:
    outname = re.sub('[^a-zA-Z0-9\n\>]', '_', inname)
    while '__' in outname:
        outname = re.sub('__', '_', outname)
    while outname[-1] == '_':
        outname = outname[:-1]
    return outname

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('infile', help="the file to clean")
  parser.add_argument('-n', action='store_true', help="whether or not to \
                          name otherwise nameless internal nodes")
  parser.add_argument('--num_ok',
                      action='store_true',
                      help="allow numbers for names")
  args = parser.parse_args()
  lines = parse_file(args.infile)
  clean = clean_seqs(lines, args.n, num_ok=args.num_ok)
  if isinstance(clean, list):
    for line in clean:
      print(line)
  else:
    print(clean)
  #print(clean)
