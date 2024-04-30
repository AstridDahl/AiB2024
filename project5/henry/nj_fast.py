import numpy as np
import time
import sys
from functools import lru_cache

time1 = time.time()

def read_control_file(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    # get the gap cost
    gap_cost = float(lines[0].strip())
    # get the alphabet
    taxa = []
    # get the scoring matrix
    dist_mat = []
    for line in lines[1:]:
        line = line.strip().split()
        taxa.append(line[0])
        dist_mat.append([float(x) for x in line[1:]])
    return gap_cost, taxa, np.array(dist_mat)

# define a basic tree class to store the result of the algorithm
class Node:
    def __init__(self, key):
        self.key = key
        self.children = []
        self.parent = None
        self.dist_to_parent = None
    def add_child(self, child):
        self.children.append(child)
    def remove_children(self, children):
        self.children = [child for child in self.children if child not in children]
    def __repr__(self):
        return str(self.key)
    
#### Helper functions #####
# calculate_r: function to calculate the r value for a given node
@lru_cache(maxsize=None)
def calculate_r(i, D):
    return np.sum(D[i, :]) / (D.shape[0] - 2)

# create_N: function to create the N matrix from the D matrix
def create_N(D):
    N = np.zeros(D.shape)
    for i in range(D.shape[0]):
        N[i, :] = np.where(i != np.arange(D.shape[0]), D[i, :] - calculate_r(i, D) - calculate_r(np.arange(D.shape[0]), D), 0)
    return N

# update_D: function to update the D matrix given the minimum indices which will be removed
def update_D(D, min_idx):
    not_min_idx = np.array([x for x in range(D.shape[0]) if x not in min_idx])
    new_distances = 0.5 * (D[min_idx[0], not_min_idx] + D[min_idx[1], not_min_idx] - D[min_idx[0], min_idx[1]])
    D_new = D[np.ix_(not_min_idx, not_min_idx)]
    D_new = np.block([[D_new, new_distances.reshape(-1, 1)], [new_distances, 0]])
    return D_new

# tree_to_newick: takes the tree output from nj and converts it to a newick format string representation
def tree_to_newick(node):
    if len(node.children) == 0:
        return f"'{node.key}'"
    else:
        return f"({','.join([tree_to_newick(x) + f':{abs(x.dist_to_parent):.3f}' for x in node.children])})" 
        
##### nj function ####
def nj(path):
    gap_cost, S, D = read_control_file(path)
    tree = Node('pseudo_root')
    for s in S:
        node = Node(s)
        node.parent = tree
        tree.add_child(node)
    while len(S) > 3:
        N = create_N(D)
        stoi = {s:i for i,s in enumerate(S)}
        itos = {i:s for s,i in stoi.items()}
        min_idx = np.unravel_index(np.argmin(N + np.eye(N.shape[0])*np.max(N)), N.shape)
        min_taxa = [itos[x] for x in min_idx]
        new_node = Node(f"{min_taxa[0]}{min_taxa[1]}")
        child1 = [x for x in tree.children if x.key == min_taxa[0]][0]
        child2 = [x for x in tree.children if x.key == min_taxa[1]][0]
        child1.parent = new_node
        child2.parent = new_node
        child1.dist_to_parent = 0.5*(D[stoi[min_taxa[0]], stoi[min_taxa[1]]] \
                                    + calculate_r(stoi[min_taxa[0]], D) \
                                    - calculate_r(stoi[min_taxa[1]], D))
        child2.dist_to_parent = D[stoi[min_taxa[0]], stoi[min_taxa[1]]] \
                                - child1.dist_to_parent
        new_node.add_child(child1)
        new_node.add_child(child2)
        tree.add_child(new_node)
        tree.remove_children([x for x in tree.children if x.key in min_taxa])
        D = update_D(D, min_idx)
        S.append(new_node.key)
        S = [x for x in S if x not in min_taxa]

    stoi = {s:i for i,s in enumerate(S)}
    itos = {i:s for s,i in stoi.items()}
    tree.children[0].dist_to_parent = 0.5*(D[stoi[tree.children[0].key], stoi[tree.children[1].key]] \
                                        + D[stoi[tree.children[0].key], stoi[tree.children[2].key]] \
                                            - D[stoi[tree.children[1].key], stoi[tree.children[2].key]])
    tree.children[1].dist_to_parent = 0.5*(D[stoi[tree.children[0].key], stoi[tree.children[1].key]] \
                                        + D[stoi[tree.children[1].key], stoi[tree.children[2].key]] \
                                            - D[stoi[tree.children[0].key], stoi[tree.children[2].key]])
    tree.children[2].dist_to_parent = 0.5*(D[stoi[tree.children[0].key], stoi[tree.children[2].key]] \
                                        + D[stoi[tree.children[1].key], stoi[tree.children[2].key]] \
                                            - D[stoi[tree.children[0].key], stoi[tree.children[1].key]])
    return tree_to_newick(tree) + ';'

def main():
    print(nj(sys.argv[1]))

if __name__ == "__main__":
    main()
    time2 = time.time()
    print(time2 - time1)
