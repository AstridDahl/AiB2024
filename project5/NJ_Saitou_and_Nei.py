# 15-04-24 
# Evt: ret floats så kun 2 decimaler. Ret parser, så tal i D ikke er strings men floats. Ret parser, så output er en liste af lister frem for en dictionary. 
# Even though the tree is unrooted, the newick representation should be rooted. An arbitrary node is selected as root node. Here, Innernode2 is used as root node. 
 
### Code to make tree in Newick format
from __future__ import annotations
from dataclasses import dataclass
from typing import Union, cast
import re
import numpy as np

def tokenize(tree: str) -> list[str]:
    """
    Extract the tokens from the text representation of a tree.

    >>> tokenize("A")
    ['A']
    >>> tokenize("(A, (B, C))")
    ['(', 'A', '(', 'B', 'C', ')', ')']
    """
    return re.findall(r'[()]|\w+', tree)


@dataclass(repr=False)
class Leaf:
    """
    A leaf in a tree.

    This will just be a string for our application.
    """

    name: str
    weight: float

    def __str__(self) -> str:
        """Simplified text representation."""
        return f"{self.name}:{self.weight}"
    __repr__ = __str__


@dataclass(repr=False)
class Node:
    """An inner node."""

    children: list[Tree]
    name: str
    weight: float

    def __str__(self) -> str:
        """Simplified text representation."""
        if type(self.weight) == float:
            return f"({','.join(str(child) for child in self.children)}){self.name}:{self.weight}"
        else:
            return f"({','.join(str(child) for child in self.children)}){self.name};"
    __repr__ = __str__


# A tree is either a leaf or an inner node with sub-trees
Tree = Union[Leaf, Node]

class EmptyStack(Exception):
    pass


class Stack(object):
    """
    Underlying data-structure is a python list.
    """
    def __init__(self):
        self.stack = []
    
    def push(self, element):
        self.stack.append(element)
    
    def get_top_element(self):
        if len(self.stack) == 0: # could also use try-except block here as on p. 552.
            raise EmptyStack()
        return self.stack[-1]
    
    def pop(self):
        if len(self.stack) == 0:
            raise EmptyStack()
        return self.stack.pop()

    def empty(self):
        return len(self.stack) == 0

    def __bool__(self): 
        return not self.empty


def parse(tree: str) -> Tree:
    """
    Parse a string into a tree.

    >>> parse("(A, (B, C))")
    (A,(B,C))
    """
    stack = Stack()
    tokens = tokenize(tree) # list of strings. 
    for token in tokens:
        if token == ')':
            # pop until ( and make (sub)tree (by making a class call to
            # Node() and give the leaf objects in a list as argument) 
            # and push the (sub)tree back onto the top of the stack.
            stack.push(token) 
            leafs =[]
            while True:
                x = stack.pop()
                if x == '(':
                    break
                elif x != ')': # if x is a leaf-object or a subtree
                    # (Node) created in line 212. 
                    leafs.append(x)
            leafs.reverse() # reverse() method reverses the list but 
            # returns None.
            subtree = Node(leafs)
            stack.push(subtree)

        elif token == '(':
            # push onto top of the stack.
            stack.push(token)

        else: # token is a string representation of a leaf.
            # create Leaf object. Push onto the top of the stack.
            stack.push(Leaf(token))
    return stack.get_top_element()

def parse_phylip(file):
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
D = parse_phylip('example_slide4.phy')
# print(D)

# Initialize. 
def initialize_S(D):
    S = []
    for key in D:
        S.append(key)
    return(S)

S = initialize_S(D)
# print(S)

def initialize_T(S):
    # Every taxon in S is a leaf in T. 
    T = []
    for i in range(len(S)):
        T.append(Leaf(name = S[i], weight = 0.00))
    return(T)

T = initialize_T(S)
# print(T)

def compute_r(i, S, D):
    '''Compute the average distance from taxon i to all other taxons in T.'''
    sum = 0
    for m in range(len(D)): # 0, 1, 2, 3, 4. 
            sum += float(D[S[i]][S[m]])  
    r = sum/(len(S) - 2)
    return(r)

r_B = compute_r(1, S, D) # r_B

def compute_n(i, j, S, D):
    '''Compute n_ij of matrix N.'''
    r_i = compute_r(i, S, D)
    r_j = compute_r(j, S, D)
    n_ij = float(D[S[i]][S[j]]) - (r_i + r_j) 
    return(n_ij)

n_BD = compute_n(1,3, S, D) # n_BD

def compute_N(S, D):
    N = {}
    for i in range(len(S)):
        N[S[i]] = {}
        for j in range(len(S)):
                if i == j:
                    N[S[i]][S[j]] = None
                else:
                    N[S[i]][S[j]] = compute_n(i, j, S, D) # Måske blot lave som liste af lister så nemmere indeksering.
    return(N)

N = compute_N(S, D)

            
def select_i(S, N): 
    n_ij = np.inf
    i = 0
    for k in range(len(S)):
        for l in range(len(S)):
            if k != l:
                if N[S[k]][S[l]] < n_ij:
                    n_ij = N[S[k]][S[l]]
                    i = k
    return(i)

#print(select_i(S, N))

def select_j(S, N):
    n_ij = np.inf
    j = 0
    for k in range(len(S)):
        for l in range(len(S)):
            if k != l:
                if N[S[k]][S[l]] < n_ij:
                    n_ij = N[S[k]][S[l]]
                    j = l
    return(j)  

#print(select_j(S, N))

def compute_weight(i, j, S, D): 
    r_i = compute_r(i, S, D)
    r_j = compute_r(j, S, D)
    w = 1/2 * (float(D[S[i]][S[j]]) + r_i - r_j)
    print(w)
    return(w)

w_ki = compute_weight(1, 3, S, D)
# print(w_ki)
w_kj = compute_weight(3, 1, S, D)
# print(w_kj)

def update_T(i, j, k, S, D, T):
    '''Add new node and new weighted edges to T.
    Remove nodes i and j from T.'''
    # Add edge between k and i to T.
    T[i].weight = compute_weight(i, j, S, D)
    # print(T[i].weight) # -0.03375
    # Add edge between k and j to T. 
    T[j].weight = compute_weight(j, i, S, D) # pass i and j in different order to get weight between k and j instead of k and i.
    # Add new node k to T. 
    T.append(Node(name = f'InnerNode{k}', weight = 0.00, children = [T[i], T[j]])) 
    # Remove child nodes from T. 
    del T[i]
    del T[j-1] # index changed. 
    return(T)

# print(update_T(1,3,0,D,T))
# S = ['A', 'B', 'C', 'D', 'E']

def update_D(i, j, k, D, S): 
    '''Remove columns and rows for (artificial) taxa i and j and add column and row for (artificial) taxa k.''' 
    D[k] = {}
    for m in range(len(S)):
        if m != i and m != j: # k not added to S yet. 
            dist = 1/2*(float(D[S[i]][S[m]]) + float(D[S[j]][S[m]]) - float(D[S[i]][S[j]])) # Fix. Only two decimal numbers. But without changing to string. 
            D[k][S[m]] = str(dist) # add row. 
            D[S[m]][k] = str(dist) # add column.
    D[k][k] = '0.00'  
    # remove rows for i and j.
    del D[S[i]]
    del D[S[j]]
    # remove columns for i and j.
    for m in range(len(S)): 
        if m != i and m != j:
            del D[S[m]][S[i]]
            del D[S[m]][S[j]]
    return(D)

#k = 0
#print(update_D(1, 3, f"InnerNode{k}", D, S))
# D = {'A': {'A': '0.00', 'C': '0.16', 'E': '0.17', 'InnerNode0': '0.13'}, 
#      'C': {'A': '0.16', 'C': '0.00', 'E': '0.11', 'InnerNode0': '0.13'}, 
#      'E': {'A': '0.17', 'C': '0.11', 'E': '0.00', 'InnerNode0': '0.13999999999999996'}, 
#      'InnerNode0': {'A': '0.13', 'C': '0.13', 'E': '0.13999999999999996', 'InnerNode0': '0.00'}}
    
def update_S(i, j, k, S):
    '''Remove (artificial) taxa i and j from S and add k.'''
    S.remove(S[i]) # remove the first matching element in S. 
    S.remove(S[j-1])
    S.append(f"InnerNode{k}")
    return(S)

# print(update_S(1,3,0,S))

def terminate(k, S, D, T): # i = 0, j = 1, m = 2. Ingen grund til at give i, j og m argumenter til funktion. 
    # Tree taxa at end instead of two, as unrooted binary tree is made. 
    # Calculate weights. 
    v_i = (float(D[S[0]][S[1]]) + float(D[S[0]][S[2]]) - float(D[S[1]][S[2]]))/2
    v_j = (float(D[S[0]][S[1]]) + float(D[S[1]][S[2]]) - float(D[S[0]][S[2]]))/2
    v_m = (float(D[S[0]][S[2]]) + float(D[S[1]][S[2]]) - float(D[S[0]][S[1]]))/2
    # Add edges to T. 
    T[0].weight = v_i
    T[1].weight = v_j
    T[2].weight = v_m
    # Add new node v to T.
    T.append(Node(name = f"InnerNode{k}", weight = '', children = [T[0], T[1], T[2]]))  
    # Delete old nodes. 
    del T[0] # idx of the next elements in the list changes when the first element is deleted. 
    del T[0] 
    del T[0]
    return(T[0]) # return the only node left in the list T. 


def NJ_saitou_and_nei(file):
    '''Takes file with distance matrix in phylip format as input and returns a 
    unrooted binary phylogenetic tree. The hierarichal clustering algorithm
    used to construct the tree is Saitou and Nei's neighbor-joining algorithm,
    which is an agglomerative hierarchical clustering algorithm. I.e., it 
    proceeds bottom-up.'''
    D = parse_phylip(file)
    # Initialization.
    S = initialize_S(D)
    T = initialize_T(S)
    # Main.
    k = 0
    while len(S) > 3:
        # Compute matrix N.
        N = compute_N(S, D)
        # Select i and j in S such that n_ij is the minimum entry in N.
        i = select_i(S, N)
        j = select_j(S, N)
        # Add new node (parent node to nodes i  and j) and weights (weight(k,i) and weight(k,j)) to T.
        # Remove old nodes from T.
        T = update_T(i, j, k, S, D, T)
        print(T)
        # Update D and S.
        D = update_D(i, j, f"InnerNode{k}", D, S)
        S = update_S(i, j, k, S)
        k += 1
    # Termination.
    T = terminate(k, S, D, T) 
    return(T)

# print(NJ_saitou_and_nei('example_slide4.phy'))

# T, D and S when termination is to be carried out: 
T = ['A:0.0', '(B:0.09999999999999998,D:0.07000000000000003)InnerNode0:0.0', '(C:0.05,E:0.06)InnerNode1:0.0']
D = {'A': {'A': '0.00', 'InnerNode0': '0.13', 'InnerNode1': '0.11000000000000001'}, 'InnerNode0': {'A': '0.13', 'InnerNode0': '0.00', 'InnerNode1': '0.07999999999999999'}, 'InnerNode1': {'A': '0.11000000000000001', 'InnerNode0': '0.07999999999999999', 'InnerNode1': '0.00'}}
S = ['A', 'InnerNode0', 'InnerNode1']

# D, T, S and N when iteration where k = 1 is to be carried out:
i = 1
j = 2
k = 0
D = {'A': {'A': '0.00', 'C': '0.16', 'E': '0.17', 'InnerNode0': '0.13'}, 'C': {'A': '0.16', 'C': '0.00', 'E': '0.11', 'InnerNode0': '0.13'}, 'E': {'A': '0.17', 'C': '0.11', 'E': 
'0.00', 'InnerNode0': '0.13999999999999996'}, 'InnerNode0': {'A': '0.13', 'C': '0.13', 'E': '0.13999999999999996', 'InnerNode0': '0.00'}}
T = ['A:0.0', 'C:0.0', 'E:0.0', '(B:0.09999999999999998,D:0.07000000000000003)InnerNode0:0.0']
S = ['A', 'C', 'E', 'InnerNode0']
N = {'A': {'A': None, 'C': -0.27, 'E': -0.27, 'InnerNode0': -0.3}, 'C': {'A': -0.27, 'C': None, 'E': -0.30000000000000004, 'InnerNode0': -0.27}, 'E': {'A': -0.27, 'C': -0.30000000000000004, 'E': None, 'InnerNode0': -0.27}, 'InnerNode0': {'A': -0.3, 'C': -0.27, 'E': -0.27, 'InnerNode0': None}}

print(NJ_saitou_and_nei('example_slide4_modified.phy'))
# (C:0.079375,(B:0.08083333333333334,D:0.08916666666666664)InnerNode1:0.050625,(A:0.044375,(E:-0.03375,F:0.13374999999999998)InnerNode0:0.09062500000000001)InnerNode2:0.03562499999999999)InnerNode3;
# Ret negativ værdi ved E. Fejl ved første iteration. Ellers korrekt. 


