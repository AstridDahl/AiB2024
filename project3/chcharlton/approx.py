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
    # make the sequences lowercase
    for header, seq in seqs.items():
        seqs[header] = seq.lower()
    return seqs

########### define function to read the Phylip-like control file ####
# this will read in the score matrix, gap cost, and alphabet
# the file will have the following format:
#   4  #gap cost
#   A  0  5  2  5 #alphabet[0], scoring matrix[0,:]
#   C  5  0  5  2 #alphabet[1], scoring matrix[1,:]
#   G  2  5  0  5 #alphabet[2], scoring matrix[2,:]
#   T  5  2  5  0 #alphabet[3], scoring matrix[3,:]

def read_control_file(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    # get the gap cost
    gap_cost = int(lines[0].strip())
    # get the alphabet
    alphabet = []
    # get the scoring matrix
    scoring_matrix = []
    for line in lines[1:]:
        line = line.strip().split()
        alphabet.append(line[0])
        scoring_matrix.append([int(x) for x in line[1:]])
    # make the alphabet lowercase
    alphabet = [x.lower() for x in alphabet]
    return gap_cost, alphabet, np.array(scoring_matrix)



######## Cost function between two nucleotides #######################
def cost(nuc1, nuc2, scoring_matrix, alphabet):
#   will use this for indexing, so it's a list with the order of alphabet
    nucleotides = alphabet
    # nucleotides = ['a','c','g','t']
#   set the scoring matrix the scoring matrix
    scoring_matrix = scoring_matrix
#   index into the scoring matrix based on the given nucleotides to get cost
    return scoring_matrix[nucleotides.index(nuc1), nucleotides.index(nuc2)]

######## define the optimality table function  #######################
def optimal(A,B, scoring_matrix, gap_cost, alphabet):
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
        T[:,0] = [0 + gap_cost*i for i in range(n+1)]
    else:
        T[:,0] = np.zeros(n+1)

    # set the 0th row of T to be 0, -gapcost, -2*gapcost, ... -n*gapcost
    # make it work for gapcost = 0
    if gap_cost != 0:
        T[0,:] = [0 + gap_cost*i for i in range(m+1)]
    else:
        T[0,:] = np.zeros(m+1)

    # fill out the table by choosing the maximum of the three options,
    # a match/mismatch defined by cost(diagonal), 
    # a gap in A (up) defined by gap cost (-1)
    # or a gap in B (left) defined by gap cost (-1)
    for i in range(1,n+1):
        for j in range(1,m+1):
            T[i,j] = min(
                T[i-1,j-1] + cost(A[i-1],B[j-1], scoring_matrix, alphabet), 
                T[i-1,j] + gap_cost, 
                T[i,j-1] + gap_cost 
            )
    return T

######## define the backtracking function  #######################
def backtrack_one(A,B, scoring_matrix, gap_cost, alphabet):
    align_A = ''
    align_B = ''
    i = len(A)
    j = len(B)
    T = optimal(A,B, scoring_matrix, gap_cost, alphabet)
    # keep iterating while we haven't reached the end of either sequence
    while i > 0 and j > 0:
        # if the score in T came from a match/mismatch...
        if T[i,j] == T[i-1,j-1] + cost(A[i-1],B[j-1], scoring_matrix, alphabet):
            align_A = A[i-1] + align_A
            align_B = B[j-1] + align_B
            i -= 1
            j -= 1
        # if the score in T came from a gap in A...
        elif T[i,j] == T[i-1,j] + gap_cost:
            align_A = A[i-1] + align_A
            align_B = '-' + align_B
            i -= 1
        # if the score in T came from a gap in B...
        else:
            align_A = '-' + align_A
            align_B = B[j-1] + align_B
            j -= 1
    return (align_A, align_B, T)

########## function to get the center string ##################################
# get optimal costs between all pairs of sequences, return the sequence that is on average closest to the others
def get_center_string(seqs):
    center_string = None
    # get the optimal costs between all pairs of sequences (not including self)
    # store the seqence key with the lowest average cost
    lowest_avg_cost = np.inf
    for key1 in seqs:
        total_cost = 0
        for key2 in seqs:
            if key1 != key2:
                A = seqs[key1]
                B = seqs[key2]
                gap_cost, alphabet, scoring_matrix = read_control_file('/Users/chcharlton/Aarhus/Spring24/AiB/AiB2024/project3/chcharlton/control.txt')
                T = optimal(A,B, scoring_matrix, gap_cost, alphabet)
                total_cost += T[-1,-1]
        avg_cost = total_cost / (len(seqs)-1)
        if avg_cost < lowest_avg_cost:
            lowest_avg_cost = avg_cost
            center_string = key1
    return center_string


################## function to generate all pairwise alignments #################
# generate alignemtns between a center string and all other sequences except itself
def generate_alignments(seqs, center_string):
    alignments = {}
    gap_cost, alphabet, scoring_matrix = read_control_file('/Users/chcharlton/Aarhus/Spring24/AiB/AiB2024/project3/chcharlton/control.txt')
    for key in seqs:
        if key != center_string:
            A = seqs[center_string]
            B = seqs[key]
            align_a, align_b, T = backtrack_one(A, B, scoring_matrix, gap_cost, alphabet)
            # store the alignments in a dictionary where each value is a numpy array of 
            # the split up strings of the alignments
            alignments[key] = np.array([list(align_a), list(align_b)])
    return alignments

################## function to extend the MSA ############################
# Takes a multiple alignment (M) and a pairwise alignemnt (A) and returns a new multiple alignment (MA)
def extend_msa(M, A):
    MA = np.empty((np.shape(M)[0]+1,0), str)
    i, j = 0, 0
    while i < np.shape(M)[1] and j < np.shape(A)[1]:
        if M[0][i] != '-' and A[0][j] != '-':
            MA = np.append(MA, np.append(M[:,i], A[-1,j])[:,None], axis=1)
            # increment both i and j because we are consuming information from both M and A
            i += 1
            j += 1
        elif M[0][i] == '-' and A[0][j] != '-':
            MA = np.append(MA, np.append(M[:,i], '-')[:, None], axis=1)
            # increment i but not j because we are consuming no information from A
            i += 1
        elif M[0][i] != '-' and A[0][j] == '-':
            MA = np.append(MA, np.append(np.array(['-']*np.shape(M)[0]), A[-1,j])[:, None], axis=1)
            # increment j but not i because we are consuming no information from M
            j += 1
        elif M[0][i] == '-' and A[0][j] == '-':
            MA = np.append(MA, np.append(M[:,i], '-')[:, None], axis=1)
            # increment i but not j because we are consuming no information from A
            i += 1
    if i < np.shape(M)[1]:
        MA = np.append(MA, np.append(M[:,i:], np.array(['-']*np.shape(A)[0])[:, None], axis=1), axis=1)
    if j < np.shape(A)[1]:
        MA = np.append(MA, np.append(np.array(['-']*np.shape(M)[0]), A[-1,j:], axis=1), axis=1)
    return MA

################## function to iteratively complete MSA ############################
# takes a dictionary of alignments and returns a multiple alignment
def get_msa(alignments):
    # get the first alignment
    M = alignments[list(alignments.keys())[0]]
    # iterate over the rest of the alignments and extend the multiple alignment
    for key in list(alignments.keys())[1:]:
        A = alignments[key]
        M = extend_msa(M, A)
    return M

# main function
def main():
    # get the sequences
    gap_cost, alphabet, scoring_matrix = read_control_file(sys.argv[1])
    seqs = read_fasta(sys.argv[2])
    # check if the sequences contain only the alphabet
    for seq in seqs.values():
        for nuc in seq:
            if nuc not in alphabet:
                raise ValueError('The sequence contains a nucleotide not in the alphabet')
    # get the sequences
    A = list(seqs.values())[0]
    B = list(seqs.values())[1]
    # get the alignment
    alignment = backtrack_one(A,B, 
                              scoring_matrix = scoring_matrix, 
                              gap_cost = gap_cost, 
                              alphabet = alphabet)
    print(alignment[2][-1][-1])
    # store the alignment in a fasta file
    if ((user_input == 'y') | (user_input == 'Y')):
        # file name: input file before its suffix + _alignment.fasta
        input_file = sys.argv[2]
        output_file = input_file[:input_file.rfind('.')] + '_alignment.fasta'
        with open(output_file, 'w') as f:
            f.write('>seq1\n' + alignment[0] + '\n')
            f.write('>seq2\n' + alignment[1] + '\n')
        
# run the main function
if __name__ == "__main__":
    main()