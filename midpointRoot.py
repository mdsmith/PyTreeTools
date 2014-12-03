#! /usr/bin/env python3

import argparse
import parsers

def midpoint_root(tree):
  new_tree = parsers.Tree(tree.tree_to_newick())
  distance_matrix = pairwise_distances(new_tree)
  max = 0
  max_distance = max([max(a) for a in distance_matrix])
  print(max_distance)
  # XXX Find the pair with the maximum distance between them, find the branch
  # spanning half the distance between them, reroot on that branch with the
  # appropriate lengths
  return new_tree

def tip_to_root(node):
  ancestry = []
  cur_node = node
  while cur_node:
    ancestry.append(cur_node)
    cur_node = cur_node.parent
  return ancestry

# Potential recursive path
def distance_between(node1, node2):
  node1_traversal = tip_to_root(node1)
  node2_traversal = tip_to_root(node2)
  if node1 == node2:
    return 0
  elif node2 in node1_traversal:
    distance = 0
    cur_node = node1
    while cur_node != node2:
      distance += float(cur_node.length)
      cur_node = cur_node.parent
    return distance
  elif node1 in node2_traversal:
    distance = 0
    cur_node = node2
    while cur_node != node1:
      distance += float(cur_node.length)
      cur_node = cur_node.parent
    return distance
  else:
    zip_rev_anc = zip(node1_traversal[::-1], node2_traversal[::-1])
    com_anc = [ancestor1 for ancestor1, ancestor2 in zip_rev_anc
                        if ancestor1 == ancestor2]
    lca = com_anc[-1]
    distance = distance_between(node1, lca) + distance_between(node2, lca)
    return distance

# Recursive:
def post_order_accumulator(node, acc, leaves_only=False):
  for n in node.children:
    post_order_accumulator(n, acc, leaves_only)
  if leaves_only:
    if len(node.children) == 0:
      acc.append(node)
  else:
    acc.append(node)

def get_leaves(tree):
  leaves = []
  post_order_accumulator(tree.root, leaves, leaves_only=True)
  return leaves

def pairwise_distances(tree):
  leaves = get_leaves(tree)
  distances = [[distance_between(i, j) for i in leaves] for j in leaves]
  return distances

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('newick_file', help="the newick file with ages and \
                                          notations")
  #parser.add_argument('-n --node', help="the node to root on")
  args = parser.parse_args()

  newick = parsers.parse_newick(args.newick_file)
  tree = parsers.Tree(newick[0].strip('\n'))
  rooted_tree = midpoint_root()
  print(tree.tree_to_newick())

