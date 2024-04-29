import numpy as np

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
    # make the alphabet lowercase
    taxa = [x.lower() for x in taxa]
    return gap_cost, taxa, np.array(dist_mat)

# define a basic tree class, should have children and parents
class Node:
    def __init__(self, key):
        self.key = key
        self.children = []
        self.parent = None
        self.dist_to_parent = None
    def add_child(self, child):
        self.children.append(child)
    # method to remove list of children
    def remove_children(self, children):
        for child in children:
            self.children.remove(child)
    def __repr__(self):
        return str(self.key)
    
#### Helper functions 
# calculate_r: function to calculate the r value for a given node
def calculate_r(i, D):
    # sum the distances of the node to all the other nodes except itself
    r = 0
    # for j in range(D.shape[0]):
    #     r += D[i,j]
    return D[i,:].sum() / (D.shape[0]-2)


# create_N: function to create the N matrix from the D matrix
def create_N(D):
    N = np.zeros(D.shape)
    for i in range(D.shape[0]):
        for j in range(D.shape[1]):
            if i != j:
                N[i,j] = (D[i,j] - calculate_r(i, D) - calculate_r(j, D))
    return N
# update_D: function to update the D matrix given the minimum indices which will be removed
def update_D(D, min_idx):
    # Get the indices that are not in min_idx
    not_min_idx = np.array([x for x in range(D.shape[0]) if x not in min_idx])
    # Calculate the new distances
    new_distances = 0.5 * (D[min_idx[0], not_min_idx] + D[min_idx[1], not_min_idx] - D[min_idx[0], min_idx[1]])
    # Create the new distance matrix by selecting the rows and columns not in min_idx
    D_new = D[np.ix_(not_min_idx, not_min_idx)]
    # Add the new distances to the last row and column
    D_new = np.block([[D_new, new_distances.reshape(-1, 1)], [new_distances, 0]])
    return D_new

# tree_to_newick: takes the tree output from nj and converts it to a newick format string representation
def tree_to_newick(node):
    if len(node.children) == 0:
        return f"'{node.key}'"
    else:
        return f"({','.join([tree_to_newick(x) + f':{abs(x.dist_to_parent):.3f}' for x in node.children])})" 
    
def nj(path):
    # read in the data
    gap_cost, S, D = read_control_file(path)

    # initialize the tree with the taxa
    tree = Node('pseudo_root')
    for s in S:
        node = Node(s)
        node.parent = tree
        tree.add_child(node)
    while len(S) > 3:
        # create the N matrix
        N = create_N(D)
        # based on S, create the dictionaries to go back and forth between indices and taxa names
        stoi = {s:i for i,s in enumerate(S)}
        itos = {i:s for s,i in stoi.items()}
        # get the index of the minimum value in the N matrix, excluding the diagonal
        # min_idx = np.unravel_index(np.argmin(D + np.eye(D.shape[0])*np.max(D)), D.shape)
        min_idx = np.unravel_index(np.argmin(N + np.eye(N.shape[0])*np.max(N)), N.shape)
        # taxa names of the minimum indices according to N
        min_taxa = [itos[x] for x in min_idx]
        # TREE CREATION
        # add a new node to the tree with the joined taxa
        new_node = Node(f"{min_taxa[0]}{min_taxa[1]}")
        # define nodes for the min_taxa, which will be children of the joined node
        # they could already be nodes with children
        child1 = [x for x in tree.children if x.key == min_taxa[0]][0]
        child2 = [x for x in tree.children if x.key == min_taxa[1]][0]
        # assign the parent of the children
        child1.parent = new_node
        child2.parent = new_node
        # assign their distances to the parent (joined node)
        child1.dist_to_parent = 0.5*(D[stoi[min_taxa[0]], stoi[min_taxa[1]]] \
                                    + calculate_r(stoi[min_taxa[0]], D) \
                                    - calculate_r(stoi[min_taxa[1]], D))
        child2.dist_to_parent = D[stoi[min_taxa[0]], stoi[min_taxa[1]]] \
                                - child1.dist_to_parent
        # add the children to the new node
        new_node.add_child(child1)
        new_node.add_child(child2)
        # add the new node to the top level of the tree
        tree.add_child(new_node)
        # remove the nodes corresponding to the joined taxa
        tree.remove_children([x for x in tree.children if x.key in min_taxa])
        D = update_D(D, min_idx)
        # add the new node to the list of taxa
        S.append(new_node.key)
        S = [x for x in S if x not in min_taxa]

    # termination: assign the distances between the last three nodes and the root
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

# use the first argument as the path to the control file
import sys
def main():
    print(nj(sys.argv[1]))
if __name__ == "__main__":
    main()