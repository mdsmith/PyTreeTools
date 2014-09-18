#! /usr/bin/env python3

import argparse

def sub_lengths(newick, new_len_dict):
  #print(newick)
  #print(new_len_dict)
  for key, value in new_len_dict.items():
    #print(newick.find(key + ":"))
    #old = newick[newick.find(key + ":"):len(key)]
    #old = newick[newick.find(key + ":"):newick.find(key + ":") + len(key)]
    start = newick.find(key + ":") + len(key) + 1
    end = min([i for i in [newick[start:].find(')'),newick[start:].find(',')]
      if i > 0]) + start
    #old = newick[start:end]
    #print(key)
    #print(len(key))
    #print(start)
    #print(end)
    #print(old)
    newick = newick[:start] + value + newick[end:]
  return newick

def get_dict(fit_file, length_formula):
  fit = open(fit_file, 'r')
  lines = fit.readlines()
  new_dict = {}
  for line in lines:
    if line[:11] == "mixtureTree":
      if '=' in line:
        name, value = line.split('=')
        param = name.split('.')[-1]
        name = name.split('.')[-2]
        #value = float(value.strip(';\n'))
        value = value.strip(';\n')
        if param in length_formula:
          new_dict[name] = value
  return new_dict

def parse_newick(newick_file):
  fh = open(newick_file, 'r')
  lines = fh.readlines()
  #newick_validate(lines)
  return lines[0]

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('newick_file', help="the newick file with ages and \
                                      notaitons")
  parser.add_argument('fit_file', help="the fit file containing possible \
                                  distance values")
  parser.add_argument('length_formula',
                      nargs=argparse.REMAINDER,
                      help="the parameter to use as the length")
  args = parser.parse_args()

  new_len_dict = get_dict(args.fit_file, args.length_formula)
  newick = parse_newick(args.newick_file)
  new_newick = sub_lengths(newick, new_len_dict)
  print(new_newick)
