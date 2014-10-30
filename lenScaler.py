#! /usr/bin/env python3

import argparse
#import sys

def scale_lengths(newick, scale):
  #print(newick)
  #print(new_len_dict)
  num = newick.count(':')
  last_end = 0
  for _ in range(num):
    #print(newick.find(key + ":"))
    #old = newick[newick.find(key + ":"):len(key)]
    #old = newick[newick.find(key + ":"):newick.find(key + ":") + len(key)]
    start = newick.find(':', last_end) + 1
    #print(key)
    #print(start)
    #print(newick)
    #print(newick[start:].find(')'))
    #print(newick[start:].find(','))
    #print(newick[start:])
    end = min([i for i in [newick[start:].find(')'),newick[start:].find(',')]
      if i > 0]) + start
    #old = newick[start:end]
    #print(key)
    #print(len(key))
    last_end = end
    #print(old)
    #print(newick[start:end])

    newick = newick[:start] + str((float(newick[start:end]) * float(scale))) + newick[end:]
  return newick

def parse_newick(newick_file):
  fh = open(newick_file, 'r')
  lines = fh.readlines()
  #newick_validate(lines)
  return lines[0].strip('\n')

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('newick_file', help="the newick file with lengths")
  parser.add_argument('scale',
                      help="the value by which to scale the lengths")
  args = parser.parse_args()

  newick = parse_newick(args.newick_file)
  new_newick = scale_lengths(newick, args.scale)
  print(new_newick)
