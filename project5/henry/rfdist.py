# ------ PACKAGES ------ #
from Bio import Phylo
from io import StringIO
from ete3 import Tree
import numpy as np
import sys
#------------------------#

########################
##       AIB P4       ##
########################

'''
Objective:
This project is about comparing evolutionary trees 
constructed using the Neighbor Joining (NJ) methods on different datasets. 
The objective is to implement an efficient algorithm for computing the
RF distance between two treesand use this implementation in an experiment. 
'''

'''
General for implementation:
- Program called rfdist 
- Input takes two (unrooted) evolutionary trees in Newick format 
- Outputs the RF distance between them.
'''

#### --------------------- PROBLEM ----------------------- ####

'''
OVERVIEW:
- Construct trees by neighbour-joining (qick-tree, rapidnj)
- Find internal branches (non-trivial splits) of tree
- Convert internal branches to bit-vectors
'''
### Compute splits (Tree -> Splits)

def newick_to_splits(newick_string):

    # Parse the Newick string to create a tree object (with ete3 library)
    # !! Could get tree from Quicktree/rapidNJ and define it here instead !!
    tree = Tree(newick_string)

    # Empty set to store splits
    splits = set()

    # Traverse the tree and collect splits
    for node in tree.traverse():
        if not node.is_leaf() and not node.is_root():
            taxa = set()
            for leaf in node:
                taxa.add(leaf.name)
            splits.add(frozenset(taxa))
    return splits


### Computate distance between splits (Splits -> Distance)

def rfdist(newick01, newick02):
    T1 = newick_to_splits(newick01)
    T2 = newick_to_splits(newick02)
    shared_splits = 0
    for split_T1 in T1:
        for split_T2 in T2:
            if split_T1 == split_T2:
                shared_splits = shared_splits +1 
    number_splits_T1 = len(T1)
    number_splits_T2 = len(T2)
    RF = (number_splits_T1 + number_splits_T2) - 2 * shared_splits
    return RF


def main():
    # get user input for the paths of the trees to be compared
    # these are the first two arguments after the scripts
    tree1_path = sys.argv[1]
    tree2_path = sys.argv[2]

    # read the trees from the files
    with open(tree1_path, 'r') as file:
        tree1 = file.read()
    with open(tree2_path, 'r') as file:
        tree2 = file.read()

    # calculate the RF distance
    rf_distance = rfdist(tree1, tree2)

    # print the result
    print(rf_distance)

if __name__ == '__main__':
    main()