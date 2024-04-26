# 21-03-24  

from Bio import Phylo
import sys
import numpy as np 

t1 = Phylo.read('tree1.new', 'newick')
t2 = Phylo.read('tree2.new', 'newick')

# Step 1-5 of Day's algorithm. All steps can be done in O(n). 
def root(t, leaf: str):
    '''Root unrooted non-binary tree at specific leaf.'''
    return(t.root_with_outgroup({'name': leaf}))

root(t1, 'seq1')
#Phylo.draw_ascii(t1)
print(t1)
#root(t2, 'seq1')
#Phylo.draw_ascii(t2)


def DF_numbering_T1(t):  
    '''Number leafs of T1 by depth-first tree traversal.'''
    ORIG = {}
    DF_number = 1
    for clade in t.find_clades():
        if clade.name != None and clade.name != 'seq1': # clade.name === None for internal node and clade.name == 'seq1' for root node in this case. 
            ORIG[clade.name] = DF_number
            DF_number += 1
    return(ORIG) 

ORIG = DF_numbering_T1(t1)
# print(ORIG)

def ORIG_to_DF(ORIG): # evt. indsæt DF_numbering_T1 her og giv t som input. 
    DF = {}
    for key in ORIG:
        DF[ORIG[key]] = key
    return(DF) # T1 og T2 har samme ORIG og DF dictionaries. 

DF = ORIG_to_DF(ORIG)
#print(DF)
#{'seq2': 1, 'seq10': 2, 'seq7': 3, 'seq3': 4, 'seq9': 5, 'seq5': 6, 'seq4': 7, 'seq8': 8, 'seq6': 9}
#{1: 'seq2', 2: 'seq10', 3: 'seq7', 4: 'seq3', 5: 'seq9', 6: 'seq5', 7: 'seq4', 8: 'seq8', 9: 'seq6'}

def label_internal_T1(t, ORIG):
    '''Label internal nodes with DF interval, i.e., interval from the minimal DF-number 
    found at a leaf of the subtree of the internal node to the maximal DF-number found
      at a leaf of the subtree of the internal node.'''
    DF_intervals = []
    for clade in t.find_clades():
        if clade.name == None and clade != t.root: # if internal node. Do so root node for which clade.name is also equl to None,
            # is not given a DF interval.  
            # Iterate over subtree of the internal node. Find min and max DF number.
            subtree = clade 
            min = np.inf 
            max = 0
            for subtree_clade in subtree.find_clades():  
                if subtree_clade.name != None and subtree_clade.name != 'seq1': # if not internal node and not root node. 
                    if min > ORIG[subtree_clade.name]:
                        min = ORIG[subtree_clade.name]
                    if max < ORIG[subtree_clade.name]:
                        max = ORIG[subtree_clade.name] 
            clade.name = '[{}, {}]'.format(min, max) # The original input tree is altered. 
            DF_intervals.append([min, max])
    return(DF_intervals) # Labeling af selve træet er in-place, så træet behøver ikke returneres.

DF_intervals1 = label_internal_T1(t1, ORIG)
Phylo.draw_ascii(t1) # Hverken weights eller navne for internal nodes shown. 


def label_internal_T2(t, ORIG):
    '''Label internal nodes of T2 with DF-intervals. An internal node is only annotated with a DF-interval, 
    if all DF-numbers of the interval are found among the leafs of the subtree of the internal node.'''
    DF_intervals = []
    for clade in t.find_clades():
        if clade.name == None and clade != t.root: # if internal node.  
            # Iterate over subtree of the internal node. Find min and max DF number.
            subtree = clade 
            min = np.inf 
            max = 0
            size = 0
            for subtree_clade in subtree.find_clades(): 
                if subtree_clade.name != None and subtree_clade.name != 'seq1': # if leaf node. 
                    size += 1
                    if min > ORIG[subtree_clade.name]: 
                        min = ORIG[subtree_clade.name]
                    if max < ORIG[subtree_clade.name]:
                        max = ORIG[subtree_clade.name]
            if max - min + 1 == size:
                clade.name = '[{}, {}]'.format(min, max) 
                DF_intervals.append([min, max])
    return(DF_intervals) # Labeling af selve træet er in-place, så træet behøver ikke returneres. 

DF_intervals2 = label_internal_T2(t2, ORIG)

def DF_intervals_function(DF_intervals1, DF_intervals2): # Works. 
    '''Return list with DF-intervals of both T1 and T2.'''
    return(DF_intervals1 + DF_intervals2)

DF_intervals = DF_intervals_function(DF_intervals1, DF_intervals2)

#def splits_function(t1, t2): # Returnerer en liste med strings fremfor en liste med lister. 
#    '''Return list of DF-intervals annotated to internal nodes of rooted, non-binary trees, T1 and T2.'''
#    splits = [] 
#    for clade in t1.find_clades():
#        if clade != t1.root and clade not in t1.get_terminals(): # if internal node.
#            if clade.name != None: 
#                splits.append(clade.name)
#    for clade in t2.find_clades():
#        if clade != t2.root and clade not in t2.get_terminals(): # if internal node.
#            if clade.name != None:
#                splits.append(clade.name)
#    return(splits)

#splits = splits_function(t1, t2)

def radix_sort(keys, m): # Works. m is the highest DF-number that a leaf has been labelled with. 
    # keys is a list containing non-trivial splits. 
    '''Sorts non-trivial splits of T1 and T2 by lexicographical order.'''
    n = len(keys)
    for j in range(2-1, -1, -1): # range(d-1, -1, -1). The number of subkeys, d, is always 2 when Day's algorithm is run. j = 1, 0. 
        # Here, subkeys at idx 
        key_buckets = [[] for bucket in range(m+1)]
        for i in range(n): # n i s the total number of non-trivial splits in T1 and T2. 
            key = keys[i]
            subkey = int(key[j])
            key_buckets[subkey].append(key)
        key_results = []
        for subkey in range(m+1):
            for key in key_buckets[subkey]:
                key_results.append(key)
        keys = key_results
    return(key_results)

m = t1.count_terminals() - 1
sorted_splits = radix_sort(DF_intervals, m)


def n_shared(sorted_splits): # Works. Takes two rooted and labelled - both on leaf nodes and internal nodes - non-binary trees. 
    '''Determine the number of shared splits, by counting the number of
    DF-intervals found in both T1 and T2, as each shared DF-interval 
    represents a shared split.'''
    shared = 0
    for i in range(len(sorted_splits)-1): # if len(splits) = 8, then i = 0..6, and i+1 = 7. 
        if sorted_splits[i] == sorted_splits[i+1]:
            shared += 1
    return(shared - 1) # minus 1 since the non-trivial split between seq1 and all the other sequences is 
# included twice in sorted_splits, due to mistakes in label_internal_T1() and label_internal_T2() such that the interval [1,9] is also included.

#print(n_shared(sorted_splits)) # should be 3. 

def n_splits(t): 
    '''Count the number of non-trivial splits in a rooted, non-binary tree.''' # Behøver ikke tælle trivial splits, 
# da disse er fælles for de to træer, givet at de to træer har n identiske leafs. 
# Number of non-trivial splits equal to the number of internal nodes minus 1. 
    count = 0
    for clade in t.find_clades():
        if clade != t.root and clade not in clade.get_terminals(): # If all descendants of the clade are not leaves 
            # and if the clade is not a leaf itself. 
            count += 1 # count internal nodes. 
    return(count - 1)

#print(n_splits(t1))
#print(n_splits(t2))

def Days(t1, t2): 
    '''Compute Robinson-Foulds distance between two unrooted non-binary trees.'''
    # The original trees are modified, i.e., in-place functions are used. 
    root(t1, 'seq1') # Step 1. Root trees. 
    root(t2, 'seq1')
    ORIG = DF_numbering_T1(t1) # Step 2. Label leaves of T1. ORIG also gives leaf labels for T2 (step 3)
    DF_intervals1 = label_internal_T1(t1, ORIG) # Label internal nodes of T1 with DF-intervals. Step 4.1.  
    DF_intervals2 = label_internal_T2(t2, ORIG) # Label internal nodes of T2 with DF-intervals. Step 4.2.
    DF_intervals = DF_intervals_function(DF_intervals1, DF_intervals2) # List of all DF-intervals in T1 and T2. Represents non-trivial splits that possible are shared by T1 and T2.
    m = t1.count_terminals() - 1 # Number of leafs in T1 minus 1. Equal to the greatest DF-number assigned to leaves of T1 and T2, 
    # as all leaves except the one turned into a root (which is also counted by get_terminals()) are given a DF-number. 
    sorted_splits = radix_sort(DF_intervals, m)
    RF_dist = n_splits(t1) + n_splits(t2) - 2 * n_shared(sorted_splits) # Step 5, find splits shared by T1 and T2. 
    return(RF_dist)

t1 = Phylo.read('tree1.new', 'newick')
t2 = Phylo.read('tree2.new', 'newick')
print('Days')
print(Days(t1, t2))