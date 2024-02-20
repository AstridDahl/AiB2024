######### import stuff ##########################################
# yes, this includes numpy. deal with it.
import numpy as np
import os
import sys

######### define a helper function to read in a fasta file ############
# this should store the sequences of the file in a dictionary

def read_fasta(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    seqs = {}
    for line in lines:
        if line[0] == '>':
            header = line.strip()
            seqs[header] = ''
        else:
            seqs[header] += line.strip()
    return seqs
########### define function to read the Phylip-like control file ####
# this will read in the score matrix, gap cost, and alphabet
def read_control_file(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    # get the alphabet
    alphabet = lines[0].strip().split()
    # get the gap cost
    gap_cost = int(lines[1].strip())
    # get the scoring matrix
    scoring_matrix = np.array([list(map(int, line.strip().split())) for line in lines[2:]])
    return alphabet, gap_cost, scoring_matrix


########### default scoring matrix ################################### 
# we won't take this as an input 
scoring_matrix = np.array([[2,0,0,0],
                           [0,2,0,0],
                           [0,0,2,0],
                           [0,0,0,2]])



######## Cost function between two nucleotides #######################
def cost(nuc1, nuc2, scoring_matrix = scoring_matrix):
#   will use this for indexing
    nucleotides = ['A', 'C', 'G', 'T']
#   set the scoring matrix the scoring matrix
    scoring_matrix = scoring_matrix
#   index into the scoring matrix based on the given nucleotides to get cost
    return scoring_matrix[nucleotides.index(nuc1), nucleotides.index(nuc2)]

######## define the optimality table function  #######################
def optimal(A,B, scoring_matrix = scoring_matrix, gap_cost = 1):
    # get table dimensions
    n = len(A)
    m = len(B)

    # initialize the table
    T = np.empty((n+1,m+1))
    # set the 0th row and column to be NaN
    T[:] = np.nan

    # set the 0th column of T to be 0, -gapcost, -2*gapcost, ... -n*gapcost
    # make it work for gapcost = 0
    if gap_cost != 0:
        # make a list of length n filled with 0, -gapcost, -2*gapcost, ... -n*gapcost
        T[:,0] = [0 - gap_cost*i for i in range(n+1)]
    else:
        T[:,0] = np.zeros(n+1)

    # set the 0th row of T to be 0, -gapcost, -2*gapcost, ... -n*gapcost
    # make it work for gapcost = 0
    if gap_cost != 0:
        T[0,:] = [0 - gap_cost*i for i in range(m+1)]
    else:
        T[0,:] = np.zeros(m+1)

    # fill out the table by choosing the maximum of the three options,
    # a match/mismatch defined by cost(diagonal), 
    # a gap in A (up) defined by gap cost (-1)
    # or a gap in B (left) defined by gap cost (-1)
    for i in range(1,n+1):
        for j in range(1,m+1):
            T[i,j] = max(
                T[i-1,j-1] + cost(A[i-1],B[j-1], scoring_matrix), 
                T[i-1,j] - gap_cost, 
                T[i,j-1] - gap_cost
            )
    return T

######## define the backtracking function  #######################
def backtrack_one(A,B, scoring_matrix = scoring_matrix, gap_cost = 1):
    align_A = ''
    align_B = ''
    i = len(A)
    j = len(B)
    T = optimal(A,B, scoring_matrix= scoring_matrix, gap_cost = gap_cost)
    # keep iterating while we haven't reached the end of either sequence
    while i > 0 and j > 0:
        # if the score in T came from a match/mismatch...
        if T[i,j] == T[i-1,j-1] + cost(A[i-1],B[j-1]):
            align_A = A[i-1] + align_A
            align_B = B[j-1] + align_B
            i -= 1
            j -= 1
        # if the score in T came from a gap in A...
        elif T[i,j] == T[i-1,j] - gap_cost:
            align_A = A[i-1] + align_A
            align_B = '-' + align_B
            i -= 1
        # if the score in T came from a gap in B...
        else:
            align_A = '-' + align_A
            align_B = B[j-1] + align_B
            j -= 1
    return (align_A, align_B)


############ user input ########################################
choice = input("Request alignment of sequences? (y/n): ")

# main function
def main():
    # get the sequences
    seqs = read_fasta(sys.argv[1])
    # get the sequences
    A = list(seqs.values())[0]
    B = list(seqs.values())[1]
    # get the alignment
    alignment = backtrack_one(A,B)
    # print the alignment
    if choice == 'y'| choice == 'Y':
        #store the alignment in a file
        with open('alignment.txt', 'w') as f:
            f.write(alignment[0] + '\n' + alignment[1])
    # print the score
    print(optimal(A,B)[-1,-1])
