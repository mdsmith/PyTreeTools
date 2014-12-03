#! /usr/bin/env python3

import argparse
import parsers
from midpointRoot import midpoint_root

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('newick_file', help="the newick file with ages and \
                                          notations")
  parser.add_argument('node', help="the node to get the average age within")
  parser.add_argument('-m --midpoint',
                      dest="midpoint",
                      help="midpoint root the tree first",
                      action="store_true")
  args = parser.parse_args()

  newick = parsers.parse_newick(args.newick_file)
  tree = parsers.Tree(newick[0].strip('\n'))
  if args.midpoint:
    rooted_tree = midpoint_root(tree)
  print(tree.avg_dist(args.node))

