#! /usr/bin/env python3

import argparse
import parsers

def midpoint_root(tree):
  # duplicate
  new_tree = parsers.Tree(tree.tree_to_newick())
  #new_tree.print_tree()
  # find pairwise leaf distances
  distance_matrix, leaves = pairwise_distances(new_tree)
  # find the two most distant leaves
  max_d = 0
  max_i = -1
  max_j = -1
  for i,r in enumerate(distance_matrix):
    for j,c in enumerate(r[i:]):
      if c > max_d:
        max_d = c
        max_i = leaves[i]
        max_j = leaves[j+i]
  # find the lowest path between them. These are most distant, so this path
  # will still be the longest in the tree
  path_between = path_through_lca(max_i, max_j)
  #for node in path_between:
    #node.print_node()
  # find the appropriate branch to break
  breakpoint = max_d/2
  #print("Max distance: " + str(max_d))
  #print("Breakpoint: " + str(breakpoint))
  cur_len = 0
  break_branch = None
  for node in path_between:
    if cur_len < breakpoint and cur_len + float(node.length) > breakpoint:
      break_branch = node
      break
    else:
      cur_len += float(node.length)
  remaining_length = breakpoint - cur_len
  new_tree.reroot(break_branch.name,
                  length=(float(break_branch.length) - remaining_length))
  #new_tree.print_tree()

  return new_tree

def tip_to_root(node):
  ancestry = []
  cur_node = node
  while cur_node:
    ancestry.append(cur_node)
    cur_node = cur_node.parent
  return ancestry

def path_through_lca(node1, node2):
  node1_traversal = tip_to_root(node1)
  node2_traversal = tip_to_root(node2)
  zip_rev_anc = zip(node1_traversal[::-1], node2_traversal[::-1])
  com_anc = [ancestor1 for ancestor1, ancestor2 in zip_rev_anc
                      if ancestor1 == ancestor2]
  lca = com_anc[-1]
  n1_to_lca = node1_traversal[:node1_traversal.index(lca)+1]
  n2_to_lca = node2_traversal[:node2_traversal.index(lca)]
  path = n1_to_lca + n2_to_lca[::-1]
  return path

def lca(node1, node2):
  node1_traversal = tip_to_root(node1)
  node2_traversal = tip_to_root(node2)
  zip_rev_anc = zip(node1_traversal[::-1], node2_traversal[::-1])
  com_anc = [ancestor1 for ancestor1, ancestor2 in zip_rev_anc
                      if ancestor1 == ancestor2]
  lca = com_anc[-1]
  return lca

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
  return distances, leaves

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('newick_file', help="the newick file with ages and \
                                          notations")
  #parser.add_argument('-n --node', help="the node to root on")
  args = parser.parse_args()

  newick = parsers.parse_newick(args.newick_file)
  tree = parsers.Tree(newick[0].strip('\n'))
  rooted_tree = midpoint_root(tree)
  print(rooted_tree.tree_to_newick())

