class Node:
  def __init__(self):
    self.children = []
    self.parent = None
    self.length = 0
    self.name = ""
  def print_node(self, depth):
    for _ in range(depth):
      print("\t", end="")
    if self.parent is not None:
      print(self.parent.name + " ", end="")
    print(self.name, ':', self.length)
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
      return node.name + ":" + node.length
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

  def tree_to_newick(self):
    self.print_newick(self.root)
    print("\n")

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
      print(node.name + ":" + node.length, end="")

  def reroot(self, node_name):
    node = self.find_node(self.root, node_name)
    #node.print_node(0)
    new_root = Node()
    new_root.children.append(node)
    new_root.children.append(node.parent)
    new_root.name = "new_RT"
    new_root.length = "0"
    for c in new_root.children:
      self.reparent(c, new_root, node)
    self.root = new_root
    self.collapse_doubles(self.root)
    new_root.print_node(0)

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
    print("name: " + node.name)
    if node.parent == new_parent:
      return
    if node is old_child:
      node.parent = new_parent
      return
    elif node.parent is None:
      node.parent = new_parent
      node.children = [c for c in node.children if c is not old_child]
      for c in node.children:
        print(c.name)
      print("new parent: " + node.parent.name)
    else:
      new_children = []
      if len(node.children) != 0:
        old_parent = node.parent
        node.parent = new_parent
        new_children.append(old_parent)
        for c in node.children:
          if c is not old_child:
            new_children.append(c)
        for c in new_children:
          print(c.name)
        print("new parent: " + new_parent.name)
        for c in new_children:
          self.reparent(c, node, node)
        node.children = [c for c in node.children if c is not old_child]
        node.children.append(old_parent)
    print("")

def make_node_ret(ns, parent):
  cur = Node()
  if ',' not in ns:
    #newNode = Node()
    #newNode.name, newNode.length = ns.split(':')
    #parent.children.append(newNode)
    cur.name, cur.length = ns.split(':')
    cur.parent = parent
  else:
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
      up = ns[ns.rfind(')') + 1:]
      if up != "":
        ns = ns[:ns.rfind(')')]
        #self.make_node(up, parent)
        if ':' in up:
          cur.name, cur.length = up.split(':')
        else:
          cur.name = up
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
