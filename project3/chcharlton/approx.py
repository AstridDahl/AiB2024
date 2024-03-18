import numpy as np
import os
import sys
from alignment_helpers import read_fasta, read_control_file, cost, optimal, backtrack_one


########## function to get the center string ##################################
# get optimal costs between all pairs of sequences, return the sequence that is on average closest to the others
def get_center_string(seqs, control_file_path):
    center_string = None
    # get the optimal costs between all pairs of sequences (excl self)
    # return the sequence key with the lowest average cost
    lowest_avg_cost = np.inf
    for key1 in seqs:
        total_cost = 0
        for key2 in seqs:
            if key1 != key2:
                A = seqs[key1]
                B = seqs[key2]
                gap_cost, alphabet, scoring_matrix = read_control_file(control_file_path)
                T = optimal(A,B, scoring_matrix, gap_cost, alphabet)
                total_cost += T[-1,-1]
        avg_cost = total_cost / (len(seqs)-1)
        if avg_cost < lowest_avg_cost:
            lowest_avg_cost = avg_cost
            center_string = key1
    return center_string


################## function to generate all pairwise alignments #################
# generate alignemtns between a center string and all other sequences except itself
def generate_alignments(seqs, center_string, control_file_path):
    alignments = {}
    gap_cost, alphabet, scoring_matrix = read_control_file(control_file_path)
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
    seqs = read_fasta(sys.argv[2])
    # get the center sequence
    center_seq = get_center_string(seqs, sys.argv[1])
    # get pairwise alignments
    pairwise_alignments = generate_alignments(seqs, center_seq, sys.argv[1])
    # get the multiple alignment
    MSA = get_msa(pairwise_alignments)
    # print the multiple alignment
    print(MSA)
        
# run the main function
if __name__ == "__main__":
    main()