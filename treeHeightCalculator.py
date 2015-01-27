#! /usr/bin/env python3
#
# Take a file, extract all trees, midpoint root, calculate height, print
# heights

import argparse
from parsers import Tree
from midpointRoot import midpoint_root
from describe import print_descriptive_statistics

def get_height(tree):
  tree = midpoint_root(tree)
  return tree.avg_dist(tree.root)

def get_trees(file_handle):
  while True:
    text_tree = file_handle.readline()
    if text_tree == '':
      break
    else:
      try:
        tree = Tree(text_tree)
        yield tree
      except ValueError:
        continue

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('tree_file', help="the file containing one or more trees")
  args = parser.parse_args()

  heights = []

  with open(args.tree_file, 'r') as fh:
    for tree in get_trees(fh):
      heights.append(get_height(tree))

  print_descriptive_statistics(heights)

