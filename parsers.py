def parse_newick(newick_file):
  fh = open(newick_file, 'r')
  lines = fh.readlines()
  #newick_validate(lines)
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

def parse_nexus(nexus_file_name):
  fh = open(nexus_file_name, 'r')
  lines = fh.readlines()
  seqs = {}
  trees = {}
  inTaxLabels = False
  inMatrix = False
  taxlabels = []
  seq_list = []
  for i,line in enumerate(lines):
    if line.strip() == "":
      continue
    if line.strip() == "TAXLABELS":
      inTaxLabels = True
      continue
    if line.strip() == "MATRIX":
      inMatrix = True
      continue
    if line.strip().startswith("TREE"):
      name, tree = line.split('=')
      name = name.replace("TREE", "").strip()
      trees[name] = tree
      continue
    if inTaxLabels:
      taxlabels = line.split(' ')
      taxlabels = [t.strip('\t') for t in taxlabels]
      taxlabels = [t.strip('\'') for t in taxlabels]
      inTaxLabels = False
    if inMatrix:
      if line.strip() == "END;":
        inMatrix = False
        continue
      seq = line.strip()
      if seq.startswith('\''):
        seq = seq.split('\'')[2].strip()
      seq = seq.replace(';', '')
      seq_list.append(seq)
  seqs = dict(zip(taxlabels, seq_list))

  return seqs, trees

def emit_fasta(seqs):
  for name,seq in seqs.items():
    print(">" + name)
    print(seq)

def fasta_to_phylip(seqs):
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

def emit_phylip(seqs, num, length):
    p_seqs, num, length = fasta_to_phylip(seqs)
    print(" " + str(num) + " " + str(length))
    for line in p_seqs:
      print(line)

def phylnames(seqs):
  import string
  import math
  sorted_keys = sorted(seqs.keys())
  p_chain = {}
  for key in sorted_keys:
    if key[1:9] not in p_chain:
      p_chain[key[1:9]] = []
    p_chain[key[1:9]].append(key)

  p_seqs = {}
  alpha = list(string.ascii_lowercase) + list(string.ascii_uppercase)
  for name, key_list in p_chain.items():
    if len(key_list) > 52:
      for i,key in enumerate(key_list):
        p_seqs[name[:-1] + alpha[math.floor(i/52)] + alpha[i%52]] = seqs[key]
    else:
      for i,key in enumerate(key_list):
        p_seqs[name + alpha[i]] = seqs[key]
  return p_seqs

# file as lines
def parse_gb(gb_file):
  import re
  members = ["DEFINITION", "ACCESSION", "FEATURES", "ORIGIN"]
  json = {}
  sequence = {}
  current_member = ""
  on = False
  for line in gb_file:
    if "//" in line:
      on = False
      json[sequence["ACCESSION"]] = sequence
      sequence = {}
      current_member = ""
    else:
      #last_current_member = current_member
      tokens = line.split()
      if len(tokens) == 0:
        continue
      last_mem = current_member
      if re.match('[A-Z]', tokens[0]):
        for member in members:
          if member == tokens[0]:
            current_member = member
        if current_member == last_mem:
          # Unknown new token detected, do not parse
          on = False
          continue
        else:
          # Known new token detected, parse for new token
          on = True
          sequence[current_member] = ' '.join(tokens[1:])
      else:
        # No new token detected, parse for old token
        if on:
          sequence[current_member] = (sequence[current_member]
                                      + ' '.join(tokens))
  return json

def json_to_fasta(  seqs_json,
                    accn = True,
                    collection_date=True,
                    definition=True,
                    clean_names=True,
                    filter_out=[]):
  #import re, sys
  import re
  from nameCleaner import clean_name
  #print(filter_out, file=sys.stderr)
  fasta_dict = {}
  for accession, json_data in seqs_json.items():
    name = ""
    if accn:
      name += accession
    if definition:
      name += "_" + json_data["DEFINITION"]

    add_tipdate = True
    if add_tipdate:
      features = json_data["FEATURES"]
      date = re.search('collection_date', features)
      if collection_date and not date:
        continue
      elif date:
        skip = False
        for word in filter_out:
          if re.search(word, features) or re.search(word, name):
            #print(word, file=sys.stderr)
            #print(accession, file=sys.stderr)
            skip = True
        if skip:
          continue
        features = features[date.end():]
        date = features.split('\"')[1]
        name += "_" + str(date)

    if clean_names:
      name = clean_name(name)
    rough_seq = json_data["ORIGIN"]
    rough_seq = rough_seq.replace(' ', '')
    seq = ''.join([c for c in rough_seq if not c.isdigit()])
    fasta_dict[name] = seq
  return fasta_dict


class Node:
  def __init__(self):
    self.children = []
    self.parent = None
    self.length = 0
    self.name = ""
  def print_node(self, depth=0):
    for _ in range(depth):
      print("\t", end="")
    if self.parent is not None:
      print(self.parent.name + " ", end="")
    print(self.name, ':', str(self.length))
    #print(len(self.children))
    for c in self.children:
      c.print_node(depth+1)

class Tree:
  #rooted = False
  #root = None
  def __init__(self, newick_string):
    self.rooted = True
    #self.root = self.make_root(newick_string)
    #self.root = Node()
    #self.make_node(newick_string.strip(';'), self.root)
    self.root = make_node_ret(newick_string.strip(';\n'), None)

  def find_node(self, cur, node_n):
    if cur.name == node_n:
      return cur
    elif cur is node_n:
      return cur
    else:
      for c in cur.children:
        find = self.find_node(c, node_n)
        if find != None:
          return find

  def avg_dist(self, node_n):
    node = self.find_node(self.root, node_n)
    if len(node.children) == 0:
      return node.length
    dist = self.subclade_distance(node, 0)
    count = self.subclade_leaves(node)
    #print(dist)
    #print(count)
    #print(dist/count)
    return dist/count

  def subclade_distance(self, node, total):
    if len(node.children) == 0:
      return total
    else:
      final_total = 0
      for c in node.children:
        final_total += self.subclade_distance(c, total + float(c.length))
      return final_total

  def subclade_leaves(self, node):
    if len(node.children) == 0:
      return 1
    else:
      final_total = 0
      for c in node.children:
        final_total += self.subclade_leaves(c)
      return final_total

  def subclade_length(self, node):
    total = 0
    for c in node.children:
      total += self.subclade_length_with_branch(c)
    return total

  def subclade_length_with_branch(self, node):
    total = 0
    for c in node.children:
      total += self.subclade_length(c)
      total += float(c.length)
    return total

  def make_root(self, ns):
    root = Node()
    #if ns[0] == "(" and ns[-1] == ")":
      #pass
    #else:
      #root.len = float(ns[ns.rfind(':'):].strip(')'))
    return root

  def make_node(self, ns, parent):
    if ',' not in ns:
      #newNode = Node()
      #newNode.name, newNode.length = ns.split(':')
      #parent.children.append(newNode)
      parent.name, parent.length = ns.split(':')
      parent.length = float(parent.length)
    else:
      pc = 0
      split = -1
      for i,c in enumerate(ns):
        if c == '(':
          pc += 1
        elif c == ')':
          pc -= 1
        elif c == "," and pc == 1:
          split = i
      if split == -1:
        exit(1)
      left = ns[1:split]
      right = ns[split + 1:ns.rfind(')')]
      up = ns[ns.rfind(')') + 1:]
      if up != "":
        #self.make_node(up, parent)
        if ':' in up:
          parent.name, parent.length = up.split(':')
        else:
          parent.name = up
      new_left = Node()
      new_right = Node()
      parent.children.append(new_left)
      parent.children.append(new_right)
      self.make_node(left, new_left)
      self.make_node(right, new_right)
    print(parent.name)

  def to_newick_string(self):
    return self.node_newick(self.root)

  def node_newick(self, node):
    if len(node.children) == 0:
      return node.name + ":" + str(node.length)
    else:
      newick = "("
      for i,c in enumerate(node.children):
        newick += self.node_newick(c)
        if i < len(node.children)-1:
          newick += ','
      newick += ')'
      return newick + str(node.name) + ":" + str(node.length)

  def print_tree(self):
    self.root.print_node(0)

  def print_newick_tree(self):
    self.print_newick(self.root)
    print("\n")

  def tree_to_newick(self):
    newick = []
    self.get_newick(self.root, newick)
    return ''.join(newick)

  def get_newick(self, node, newick):
    if len(node.children) == 0:
      newick.append(node.name + ":" + str(node.length))
    else:
      newick.append('(')
      for i,c in enumerate(node.children):
        self.get_newick(c, newick)
        if i < len(node.children)-1:
          newick.append(',')
      newick.append(')')
      newick.append(node.name + ":" + str(node.length))

  def print_newick(self, node):
    if len(node.children) == 0:
      print(node.name + ":" + node.length, end="")
    else:
      print('(', end="")
      for i,c in enumerate(node.children):
        self.print_newick(c)
        if i < len(node.children)-1:
          print(',', end="")
      print(')', end="")
      print(node.name + ":" + str(node.length), end="")

  # current hypothesis, need to make a new dummy root, make a second node for
  # the node to be split, hitch them both up to the new dummy root, then
  # destroy the old dummy root
  def reroot(self, node_name, pointer=None, length=None):
    if node_name == self.root.name:
      return
    if pointer:
      node = pointer
    else:
      node = self.find_node(self.root, node_name)
    #node.print_node(0)
    if length:
      length_to_old_parent = float(node.length) - length
      node.length = length
      #node.parent.length = length_to_old_parent
    new_root = Node()
    new_root.name = "new_RT"
    new_root.length = length_to_old_parent
    new_root.children.append(node)

    #new_twig = Node()
    #new_twig.name = node.name + "_new"
    #new_twig.length = length_to_old_parent

    #new_twig.children.append(node.parent)

    #new_root.children.append(new_twig)
    #new_root.children.append(node)

    #trav = [node]

    #node.print_node()

    #self.root.print_node()

    self.replace_child(node.parent, node, new_root)
    new_root.parent = node.parent
    node.parent = new_root

    #self.root.print_node()

    trav = [new_root]
    while trav[-1].parent:
      trav.append(trav[-1].parent)

    #print(len(trav))
    #for t in trav:
      #print(t.print_node())

    new_parent = None
    new_length = 0
    for n in trav:
      new_length = self.rotate(n, new_parent, new_length)
      new_parent = n
    #print("Residual: " + new_length)
    #new_twig.parent = new_root

    #for c in new_twig.children:
      #self.reparent(c, new_twig, node)
    new_root.parent = None
    self.root = new_root
    #self.root.print_node()

    self.collapse_doubles(self.root)
    #new_root.print_node(0)

  def replace_child(self, parent, old_child, new_child):
    new_children = []
    for c in parent.children:
      if c is not old_child:
        new_children.append(c)
    new_children.append(new_child)
    parent.children = new_children

  def rotate(self, node, new_parent, length=None):
    '''
    print("Name: " + node.name)
    if node.parent != None:
      print("Parent: " + node.parent.name)
    else:
      print("Parent: None")
    print("Length: " + str(node.length))
    print("Children:")
    for c in node.children:
      print(c.name)
    print("")
    '''

    if length == None:
      length = new_parent.length
    old_length = node.length

    old_parent = node.parent

    old_children = []
    for c in node.children:
      if c is not new_parent:
        old_children.append(c)

    node.length = length
    node.parent = new_parent
    node.children = []
    if old_parent:
      node.children.append(old_parent)
    node.children += old_children

    '''
    print("Name: " + node.name)
    if node.parent != None:
      print("Parent: " + node.parent.name)
    else:
      print("Parent: None")
    print("Length: " + str(node.length))
    print("Children:")
    for c in node.children:
      print(c.name)
    print("")
    '''
    return old_length

  def collapse_doubles(self, node):
    if len(node.children) == 0:
      return
    elif len(node.children) >= 2:
      for c in node.children:
        self.collapse_doubles(c)
      return
    else:
      node.children[0].parent = node.parent
      node.parent.children = [c for c in node.parent.children if c is not
          node] + node.children
      self.collapse_doubles(node.children[0])

  def reparent(self, node, new_parent, old_child):
    #print("name: " + node.name)
    if node.parent == new_parent:
      return
    if node is old_child:
      node.parent = new_parent
      return
    elif node.parent is None:
      node.parent = new_parent
      node.children = [c for c in node.children if c is not old_child]
      #for c in node.children:
        #print(c.name)
      #print("new parent: " + node.parent.name)
    elif new_parent in node.children or old_child in node.children:
      '''
      print("Node: " + node.name)
      print("Old Parent: " + node.parent.name)
      print("Old Children: ", end="")
      for c in node.children:
        print(c.name)
      '''
      old_parent = node.parent
      old_children = []
      for c in node.children:
        if c is not old_child:
          old_children.append(c)
      old_siblings = []
      for c in node.parent.children:
        if c is not node:
          old_siblings.append(c)

      node.parent = new_parent
      node.parent.children += old_children
      node.children = []
      node.children.append(old_parent)
      node.children += old_siblings
      '''
      print("Node: " + node.name)
      print("Parent: " + node.parent.name)
      print("Children: ", end="")
      for c in node.children:
        print(c.name)
      '''
      #exit(0)
      for c in node.children:
        self.reparent(c, node, node)
    else:
      node.parent = new_parent

      '''
      if len(node.children) != 0:
        old_parent = node.parent
        node.parent = new_parent
        new_children.append(old_parent)
        for c in node.children:
          if c is not old_child:
            new_children.append(c)
        #for c in new_children:
          #print(c.name)
        #print("new parent: " + new_parent.name)
        for c in new_children:
          new_parent.children.append(c)
          self.reparent(c, node, node)
        node.children = [c for c in node.children if c is not old_child]
        node.children.append(old_parent)
    #print("")
    '''

def make_node_ret(ns, parent):
  cur = Node()
  if ',' not in ns:
    #newNode = Node()
    #newNode.name, newNode.length = ns.split(':')
    #parent.children.append(newNode)
    #print(ns)
    cur.name, cur.length = ns.split(':')
    cur.parent = parent
  else:
    #print(ns)
    pc = 0
    #split = -1
    splits = []
    for i,c in enumerate(ns):
      if c == '(':
        pc += 1
      elif c == ')':
        pc -= 1
      elif c == "," and pc == 1:
        #split = i
        splits.append(i)
    #if split == -1:
    if len(splits) == 0:
      print("Split point not found...")
      print(ns)
      exit(1)
    elif len(splits) == 1:
      split = splits[0]
      left = ns[1:split]
      right = ns[split + 1:ns.rfind(')')]
      up = ns[ns.rfind(')') + 1:]
      if up != "":
        #self.make_node(up, parent)
        if ':' in up:
          cur.name, cur.length = up.split(':')
        else:
          cur.name = up
      cur.children.append(make_node_ret(left, cur))
      cur.children.append(make_node_ret(right, cur))
    else:
      #print("polytomy")
      up = ns[ns.rfind(')') + 1:]
      if up != "":
        #ns = ns[:ns.rfind(')')]
        #self.make_node(up, parent)
        if ':' in up:
          cur.name, cur.length = up.split(':')
        else:
          cur.name = up
      ns = ns[:ns.rfind(')')]
      chunks = []
      full = ns
      last = 0
      splits.append(len(ns))
      for split in splits:
        chunks.append(full[last+1:split])
        last = split
      for chunk in chunks:
        cur.children.append(make_node_ret(chunk, cur))
    cur.parent = parent
  #print(cur.name)
  return cur
