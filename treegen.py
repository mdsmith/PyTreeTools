#! /usr/bin/env python

import sys

def get_tree(num_taxa):
    result = gen_newick([a for a in range(1, int(num_taxa) + 1)],
                        (2*int(num_taxa)) - 1)
    return  ("(" + str(result) + ")")

def gen_newick(part, number):
    if len(part) == 2:
        return  ("("
                + str(part[0])
                + ","
                + str(part[1])
                + ")"
                + str(number))
    else:
        return  ("("
                + gen_newick(   part[:int(len(part)/2)],
                                number - int(len(part)/2))
                + ","
                + gen_newick(part[int(len(part)/2):], number - 1)
                + ")"
                + str(number))

def get_unrooted_iteratively(num_taxa):
    tree = ["(", "N", ");"]
    covered = 1
    index = 0
    while covered < num_taxa:
        if tree[index] == "N" and (covered == 1 or covered == 2):
            tree[index] = "N"
            tree.insert(index + 1, ",")
            tree.insert(index + 2, "N")
            covered += 1
            index += 3
        elif tree[index] == "N":
            tree[index] = "("
            tree.insert(index + 1, "N")
            tree.insert(index + 2, ",")
            tree.insert(index + 3, "N")
            tree.insert(index + 4, ")")
            tree.insert(index + 5, "I")
            covered += 1
            index += 6
        else:
            index += 1
        index = index % len(tree)

    leafI = 1
    intI = num_taxa + 1
    for i,c in enumerate(tree):
        if c == "N":
            tree[i] = str(leafI)
            leafI += 1
        if c == "I":
            tree[i] = str(intI)
            intI += 1

    return "".join(tree)


def get_rooted_iteratively(num_taxa):
    tree = ["(", "N", ");"]
    covered = 1
    index = 0
    while covered < num_taxa:
        if tree[index] == "N":
            tree[index] = "("
            tree.insert(index + 1, "N")
            tree.insert(index + 2, ",")
            tree.insert(index + 3, "N")
            tree.insert(index + 4, ")")
            tree.insert(index + 5, "I")
            covered += 1
            index += 6
        else:
            index += 1
        index = index % len(tree)

    leafI = 1
    intI = num_taxa + 1
    for i,c in enumerate(tree):
        if c == "N":
            tree[i] = str(leafI)
            leafI += 1
        if c == "I":
            tree[i] = str(intI)
            intI += 1

    return "".join(tree)

def get_unrooted_tree(num_taxa):
    result = gen_unrooted_newick(   [a for a in range(1, int(num_taxa) + 1)],
                                    (2*int(num_taxa))-1,
                                    (2*int(num_taxa))-3)
    return("(" + result + ")")

def gen_unrooted_newick(part, number, max):
    if len(part) == 2:
        if number > max:
            return  str(part[0]) + "," + str(part[1])
        else:
            return  "(" + str(part[0]) + "," + str(part[1]) + ")" + \
                    str(number)
    else:
        if number > max:
            return  (gen_unrooted_newick(part[:int(len(part)/2)], number -
                    int(len(part)/2), max)
                    + ","
                    + gen_unrooted_newick(  part[int(len(part)/2):],
                                            number - 1,
                                            max))
        else:
            return  ("("
                    + gen_unrooted_newick(  part[:int(len(part)/2)],
                                            number - int(len(part)/2),
                                            max)
                    + ","
                    + gen_unrooted_newick(part[int(len(part)/2):], number - 1, max)
                    + ")"
                    + str(number))


if __name__ == "__main__":
    if len(sys.argv) == 2:
        numBranches = sys.argv[1]
        type = "unrooted"
    elif len(sys.argv) == 3:
        numBranches = sys.argv[1]
        type = sys.argv[2]
    else:
        print(  "Valid usage:\n\t- treegen <number of taxa> [<type>]" \
                "\n\t\t-<type>: rooted or unrooted",
                file=sys.stderr)
        exit(1)
    if type == "rooted":
        #print(get_tree(numBranches))
        print(get_rooted_iteratively(int(numBranches)))
    elif type == "unrooted":
        #print(get_unrooted_tree(numBranches))
        print(get_unrooted_iteratively(int(numBranches)))
