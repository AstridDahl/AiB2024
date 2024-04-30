# ------ PACKAGES ------ #
from Bio import Phylo
import numpy as np
import re
import io
#------------------------#

########################
##       AIB lolol    ##
########################

##--------------------------------------------------------------##
def parse_phylip(file):
    '''Takes phylip file with distance matrix and returns dictionary with distance matrix.'''
    S = [] 
    D = {}
    f = open(file, 'r')
    dist_matrix = []
    for line in f:
        lst = line.strip().split()
        if len(lst) != 1: # Not add integer in first line to S. 
            dist_matrix.append(lst[1:])
            S.append(lst[0])
    for i in range(len(S)):
        D[S[i]] = {}
        for j in range(len(S)):
            D[S[i]][S[j]] = dist_matrix[i][j]
    return(D)

##----------------------------------------------------------------##

### Make a 'None' matrix (SxS):
def none_matrix(S):
    matrix = {}
    for row_label in S:
        row = {}
        for col_label in S:
            row[col_label] = None
        matrix[row_label] = row
    return matrix

### Calculate value i 'nij':
# Get row sums:
def row_sums(distance_matrix):
    row_sums = {}
    for key, row in distance_matrix.items():
        row_sum = sum(float(value) for value in row.values())
        row_sums[key] = row_sum
    return row_sums


# Compute ri 
def r(sequence_name, distance_matrix):
    S = list(distance_matrix.keys())
    rs_d = row_sums(distance_matrix)
    r_seq = (rs_d[sequence_name])/(len(S)-2)
    return r_seq


# Compute nij 
def nij_calc(sequence_1, sequence_2, distance_matrix):
    r_1 = r(sequence_1, distance_matrix)
    r_2 = r(sequence_2, distance_matrix)
    dij = float(distance_matrix[sequence_1][sequence_2])
    nij = dij - (r_1 + r_2)
    return nij


### Make N matrix with nij values:
def N(distance_matrix):
    S = list(distance_matrix.keys())
    N = none_matrix(S)
    for i in S:
        for j in S:
            nij = nij_calc(i, j, distance_matrix)
            N[i][j] = nij
    return N


### Find minimum entry in N (that's not within the same sequence):

def min_entry_N(N_matrix):
    n = len(N_matrix)
    S = list(N_matrix.keys())
    min_value = float("inf")
    for i in S:
        for j in S:
            if i != j and N_matrix[i][j] < min_value:
                min_value = N_matrix[i][j]
    return min_value

# Find out which species are connected to min. entry

def find_keys_for_value(nm, min_val):
    keys_with_value = []

    for key_outer, inner_dict in nm.items():
        for key_inner, value in inner_dict.items():
            if value == min_val:
                keys_with_value.append((key_outer, key_inner))

    return keys_with_value[0]

### Add 'i' and 'j' together (with min. nij) in node 'ij':

# Cluster naming:
def cluster_naming(N, min_val):
    keys = find_keys_for_value(N, min_val)
    cluster = ''
    for e in keys: 
        cluster += str(e) 
    return cluster

# Calculate weigths for branches to new cluster:

def weight(i, j, distance_matrix):
    r_i = r(i, distance_matrix)
    r_j = r(j, distance_matrix)
    dij = float(distance_matrix[i][j])
    cluster_i_dist = 1/2 * (dij + r_i - r_j)
    cluster_j_dist = 1/2* (dij + r_j - r_i)
    return cluster_i_dist, cluster_j_dist


# New column and row distances in D (from k to new cluster) (k = all other sequences)
                
# D_cluster_k = 1/2(D_i_k + D_j_k - Dij)
def new_D(i, j, cluster, D, S): 
    # Add key is D with the 'cluster' name
    DD = D
    DD[cluster] = {}
    # All seqs. in S that are not 'B' or 'D' ('A','C','E').
    for m in S:
        if m != i and m != j: 
            dist = 1/2*(float(DD[i][m]) + float(DD[j][m]) - float(DD[i][j])) # 1/2(D['B']['A'] + D['D']['A'] - D['B']['D'])
            DD[cluster][m] = str(dist) # D['BD']['A'] = dist (row)
            DD[m][cluster] = str(dist) # D['A']['BD'] = dist (column)

    # Specify that's = 0.00
    DD[cluster][cluster] = '0.00'  

    # remove rows 'i' and 'j' ('A' and 'D').
    del DD[i] # D['A']
    del DD[j] # D['D']

    # remove columns for i and j.
    for m in S: 
        if m != i and m != j:
            del DD[m][i] # D['A']['B']
            del DD[m][j] # D['A']['D']
    
    return DD


# Update S
def new_S(i, j, cluster, S):
    S.remove(i) # remove the first matching element in S. 
    S.remove(j)
    S.append(cluster)
    return(S)


# Make dict for tree
def make_tree(i, j, w1, w2):
    dict = {i: w1, j: w2}
    return dict


###  ---------- Termination -------------- ##
def termination(D, S):
    v_i = (float(D[S[0]][S[1]]) + float(D[S[0]][S[2]]) - float(D[S[1]][S[2]]))/2
    v_j = (float(D[S[0]][S[1]]) + float(D[S[1]][S[2]]) - float(D[S[0]][S[2]]))/2
    v_m = (float(D[S[0]][S[2]]) + float(D[S[1]][S[2]]) - float(D[S[0]][S[1]]))/2
    dict = {S[0]: v_i, S[1]: v_j, S[2]: v_m}
    return dict

## ---------------------- Newickify a nested dictionary ---------------- ##

def newickify(node_to_children, root_node) -> str:
    visited_nodes = set()

    def newick_render_node(name, distance: float) -> str:
        assert name not in visited_nodes, "Error: The tree may not be circular!"

        if name not in node_to_children:
            # Leafs
            return F'{name}:{distance}'
        else:
            # Nodes
            visited_nodes.add(name)
            children = node_to_children[name]
            children_strings = [newick_render_node(child, children[child]) for child in children.keys()]
            children_strings = ",".join(children_strings)
            return F'({children_strings}){name}:{distance}'

    newick_string = newick_render_node(root_node, 0) + ';'

    # Ensure no entries in the dictionary are left unused.
    assert visited_nodes == set(node_to_children.keys()), "Error: some nodes aren't in the tree"

    return newick_string

### ---------- DRAW TREE FROM NEWICK --------------- ###

def newick_to_tree(newick_str):
    # Create a file-like object from the Newick string
    file_handle = io.StringIO(newick_str)
    
    # Parse the Newick string
    tree = Phylo.read(file_handle, "newick")

    # Draw the tree
    Phylo.draw(tree)

################################################################## NJ ##########################################################################

def nj(phy_file):

    # Input:
    dm = parse_phylip(phy_file)

    # List of the taxons:
    S = list(dm.keys())

    # Nested dict:
    newick = {}

    while len(S) > 3:
        ### Calculations:
        # Make N from D
        N_matrix = N(dm)

        # Fin min. entry
        min_value = min_entry_N(N_matrix)

        # Find 'i' and 'j' that corresponds to min. entry
        min_pair = find_keys_for_value(N_matrix, min_value)
        i, j = min_pair[0], min_pair[1]

        # Give a cluster name
        cluster_name = cluster_naming(N_matrix, min_value)
        
        # Calculate edges
        weights = weight(i, j, dm)
        w1, w2 = weights[0], weights[1]

        ### Tree contruction:
        # Make dict
        dict = make_tree(i, j, w1, w2)
 
        # Add dict to newick
        newick[cluster_name] = dict

        ### New D and S:
        # Calculate new D
        dm = new_D(i, j, cluster_name, dm, S)
    
        # Reduce S
        S = new_S(i, j, cluster_name, S)

    # Termination:
    print("last dm", dm)
    last_node = 'v'
    last_dd = termination(dm, S)
    newick[last_node] = last_dd

    return newick

print(nj('example_slide4.phy'))

tree = {'BD': {'B': 0.09999999999999998, 'D': 0.07000000000000003}, 
 'CE': {'C': 0.05, 'E': 0.06}, 
 'v': {'A': 0.08000000000000002, 'BD': 0.04999999999999999, 'CE': 0.03}}

print(newickify(tree, root_node = 'v'))

newick_string = "(((leaf1:12,leaf2:32)a:2,leaf3:21,leaf4:3)b:3,(leaf5:5,leaf6:7)c:5)root:0;"
newick_tree = "(A:0.08000000000000002,(B:0.09999999999999998,D:0.07000000000000003)BD:0.04999999999999999,(C:0.05,E:0.06)CE:0.03)v:0;"

print(newick_to_tree(newick_tree))