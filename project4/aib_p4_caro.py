# ------ PACKAGES ------ #
from Bio import Phylo
from io import StringIO
from ete3 import Tree
import numpy as np
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

tree_01 = "((seq2:0.18258,(seq10:0.14710,seq7:0.15261):0.02441):0.00702,seq3:0.20934,(((seq5:0.30207,seq4:0.24793):0.02587,(seq1:0.20160,(seq8:0.12843,seq6:0.12375):0.06209):0.01638):0.02352,seq9:0.17922):0.00628);"
tree_02 = "((seq1:0.18258,(seq10:0.14710,seq7:0.15261):0.02441):0.00702,seq3:0.20934,(((seq5:0.30207,seq4:0.24793):0.02587,(seq2:0.20160,(seq8:0.12843,seq6:0.12375):0.06209):0.01638):0.02352,seq9:0.17922):0.00628);"

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

# rf_test_dist = rfdist(tree_01, tree_02)  
# print("RF test distance is:", rf_test_dist, ", and it should be 8.")

#### --------------------- EXPERIMENT ----------------------- ####  

### Contruct tree (Quicktree/rapidNJ) (Newick -> Tree)

'''
Did that in the WSL.
Required linux system to run
'''

### Calculate RF distance:

## Clustal:
# quicktree:
clustal_quick_path = 'e1_newick/clustal_quicktree.newick'
with open(clustal_quick_path, 'r') as file:
    clustal_quick = file.read()

# rapidnj:
clustal_rapidnj_path = 'e1_newick/clustal_rapidnj.newick'
with open(clustal_rapidnj_path, 'r') as file:
    clustal_rapidnj = file.read()

## Muscle:


## Kalign: 

# RF-distances:
rf_clustal_rapid_clustal_quick = rfdist(clustal_quick, clustal_rapidnj) 

print(rf_clustal_rapid_clustal_quick)






