
# EVAL questions (in the end of the file)

# For outputting nice looking matrices:
def print_matrix(matrix):
    for row in matrix:
        print(" ".join(map(str, row)))
#-----------------------------------------#
        
#----------------#
#                #
# GLOBAL LINEAR  #
#                #
#----------------#

def needleman_wunsch(seq1, seq2, score_matrix, gap_penalty):
    # Initialize the scoring matrix
    n = len(seq1)
    m = len(seq2)
    S = [[0] * (m + 1) for _ in range(n + 1)]

    # Initialize the first row and first column of the matrix
    for i in range(1, n + 1):
        S[i][0] = S[i - 1][0] + gap_penalty
    for j in range(1, m + 1):
        S[0][j] = S[0][j - 1] + gap_penalty

    # Fill the scoring matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # Calculate the scores for the three possible operations
            diagonal = S[i - 1][j - 1] + score_matrix[seq1[i-1]][seq2[j-1]]
            up = S[i - 1][j] + gap_penalty
            left = S[i][j - 1] + gap_penalty

            # Choose the maximum score
            S[i][j] = min(diagonal, up, left)

    # Traceback to find the alignment
    alignment1 = ''
    alignment2 = ''
    i = n
    j = m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and S[i][j] == S[i - 1][j - 1] + score_matrix[seq1[i - 1]][seq2[j - 1]]:
            alignment1 = seq1[i - 1] + alignment1
            alignment2 = seq2[j - 1] + alignment2
            i -= 1
            j -= 1
        elif i > 0 and S[i][j] == S[i - 1][j] + gap_penalty:
            alignment1 = seq1[i - 1] + alignment1
            alignment2 = '-' + alignment2
            i -= 1
        else:
            alignment1 = '-' + alignment1
            alignment2 = seq2[j - 1] + alignment2
            j -= 1

    return alignment1, alignment2, S[n][m]


#----------------#
#                #
# GLOBAL AFFINE  #
#                #
#----------------#


####### IMPORTS ########################################
import numpy as np
import os
import sys

###### AlIGNMENT FUNCTIONS ################################
def empty_matrix(r, c): 
    m = []
    for i in range(r):
        m.append([])
        for j in range(c):
            m[i].append(None)
    return m

def prepare_matrix_D(r, c, a, b):
    pm = empty_matrix(r, c)
    pm[0][0] = 0
    for j in range(1, c): 
        pm[0][j] = a+b*(j-1) 
    for i in range(1, r):
        pm[i][0] = a+b*(i-1)
    return pm

def prepare_matrix_E(r, c):
    pm = empty_matrix(r, c)
    for i in range(c):
        pm[0][i] = np.inf 
    return pm


def prepare_matrix_F(r, c):
    pm = empty_matrix(r, c)
    for i in range(r):
        pm[i][0] = np.inf
    return pm 

substitution_matrix2 = {'A': {'A': 0, 'C': 5, 'G': 2, 'T': 5}, 
                      'C': {'A': 5, 'C': 0, 'G': 5, 'T': 2}, 
                      'G': {'A': 2, 'C': 5, 'G': 0, 'T': 5}, 
                      'T': {'A': 5, 'C': 2, 'G': 5, 'T': 0}}
def compute_cost(seq1, seq2, scores, a, b):
    r = len(seq1)+1
    c = len(seq2)+1
    D = prepare_matrix_D(r, c, a, b) # Matrix called S in slides. 
    E = prepare_matrix_E(r, c) # Matrix called D in slides. 
    F = prepare_matrix_F(r, c) # Matrix called I in slides. 
    for i in range(1, r): # row index. i from 1 to r-1, i.e., i from 1 to n1.
        for j in range(1, c): # col index. j from 1 to n2. 
            E[i][j] = min(D[i-1][j]+a, E[i-1][j]+b) # Never j-1 in E, therefore, it doesn't matter that E[i,0] is not filled out for i = 1..n1. 
            F[i][j] = min(D[i][j-1]+a, F[i][j-1]+b) # Never i-1 in F, therefore, it doesn't matter that F[0,j] is not filled out for j = 1..n2. 
            D[i][j] = min(E[i][j], F[i][j], D[i-1][j-1] + scores[seq1[i-1].capitalize()][seq2[j-1].capitalize()]) # i-1 and j-1 
            # when indexing sequences, as seq1[n1] and seq2[n2] are out of range in Python, where the first index is 0. 
    return D[r-1][c-1] # minus 1, since cell (r, c) is not found in matrix D, as the first cell in the matrix is given by (0,0) and not (1,1). 

def fill_matrices(seq1, seq2, scores, a, b): 
    r = len(seq1)+1
    c = len(seq2)+1
    D = prepare_matrix_D(r, c, a, b)
    E = prepare_matrix_E(r, c)
    F = prepare_matrix_F(r, c)
    for i in range(1, r):
        for j in range(1, c):  
            E[i][j] = min(D[i-1][j]+a, E[i-1][j]+b)
            F[i][j] = min(D[i][j-1]+a, F[i][j-1]+b)  
            D[i][j] = min(E[i][j], F[i][j], D[i-1][j-1] + scores[seq1[i-1].capitalize()][seq2[j-1].capitalize()]) 

    matrices = {'D': D, 'E': E, 'F': F}
    return  matrices 

def traceback(seq1, seq2, matrices, scores, a, b): 
  
    D = matrices['D']
    E = matrices['E']
    F = matrices['F']

    aligned1 = ''
    aligned2 = ''
    i = len(seq1)
    j = len(seq2)

    while i > 0 and j > 0: # > evalueres før 'and'. 
        # last character of seq1 is seq1[i-1]. last character of seq2 is seq2[j-1].
        # If opening gap.  
        if D[i][j] == D[i-1][j-1] + scores[seq1[i-1].capitalize()][seq2[j-1].capitalize()]:
            aligned1 = seq1[i-1] + aligned1
            aligned2 = seq2[j-1] + aligned2
            i -= 1
            j -= 1
        elif D[i][j] == D[i][j-1] + a: # left. a is opening cost. 
            aligned1 = '-' + aligned1
            aligned2 = seq2[j-1] + aligned2
            j -= 1
        elif D[i][j] == D[i-1][j] + a: # right. a is opening cost. 
            aligned1 = seq1[i-1] + aligned1
            aligned2 = '-' + aligned2
            i -= 1   
        # If extension gap. 
        else: 
            k = 2
            while k <= i and k <= j: # '<=' evalueres før 'and'.
                if D[i][j] == D[i][j-k] + a + b*(k-1):
                    aligned1 = '-'*k + aligned1
                    aligned2 = seq2[j-k:j] + aligned2 
                    j -= k
                    break
                elif D[i][j] == D[i-k][j] + a + b*(k-1):
                    aligned1 = seq1[i-k:i] + aligned1 
                    aligned2 = '-'*k + aligned2
                    i -= k
                    break
                k += 1 
            while k <= j:
                if D[i][j] == D[i][j-k] + a + b*(k-1):
                    aligned1 = '-'*k + aligned1
                    aligned2 = seq2[j-k:j] + aligned2 
                    j -= k
                    break
                k += 1
            while k <= i:
                if D[i][j] == D[i-k][j] + a + b*(k-1):
                    aligned1 = seq1[i-k:i] + aligned1 
                    aligned2 = '-'*k + aligned2 
                    i -= k
                    break
                k += 1     
    if i > 0: # j = 0. 
        aligned1 = seq1[:i] + aligned1
        aligned2 = '-'*i + aligned2
    elif j > 0: # i = 0. 
        aligned1 = '-'*j + aligned1
        aligned2 = seq2[:j] + aligned2
    return [aligned1, aligned2]

def align(seq1, seq2, scores, a, b):
    n = len(seq1)
    m = len(seq2)
    matrices = fill_matrices(seq1, seq2, scores, a, b)
    D = matrices['D']
    alignment = traceback(seq1, seq2, matrices, scores, a, b)
    return alignment, D[n][m]


#----------------#
#                #
# EVAL QUESTIONS #
#                #
#----------------#

############### Information:
d = {1 : "tatggagagaataaaagaactgagagatctaatgtcgcagtcccgcactcgcgagatactcactaagaccactgtggaccatatggccataatcaaaaag".upper(), 
     2 : "atggatgtcaatccgactctacttttcctaaaaattccagcgcaaaatgccataagcaccacattcccttatactggagatcctccatacagccatggaa".upper(), 
     3 : "tccaaaatggaagactttgtgcgacaatgcttcaatccaatgatcgtcgagcttgcggaaaaggcaatgaaagaatatggggaagatccgaaaatcgaaa".upper(),
     4 : "aaaagcaacaaaaatgaaggcaatactagtagttctgctatatacatttgcaaccgcaaatgcagacacattatgtataggttatcatgcgaacaattca".upper(),
     5 : "atgagtgacatcgaagccatggcgtctcaaggcaccaaacgatcatatgaacaaatggagactggtggggagcgccaggatgccacagaaatcagagcat".upper()}
score_matrix = {
    'A': {'A': 0, 'C': 5, 'G': 2, 'T': 5},
    'C': {'A': 5, 'C': 0, 'G': 5, 'T': 2},
    'G': {'A': 2, 'C': 5, 'G': 0, 'T': 5},
    'T': {'A': 5, 'C': 2, 'G': 5, 'T': 0}
}
gap_open_penalty = 10
gap_extension_penalty = 5
gap_penalty = 5
#############################

# Q1
print("report, Q1, global linear")
ss1 = "tatggagagaataaaagaactgagagatctaatgtcgcagtcccgcactcgcgagatactcactaagaccactgtggaccatatggccataatcaaaaag".upper()
ss2 = "atggatgtcaatccgactctacttttcctaaaaattccagcgcaaaatgccataagcaccacattcccttatactggagatcctccatacagccatggaa".upper()
print(needleman_wunsch(ss1, ss2, score_matrix, gap_penalty))

# Q2
print("report, Q2, global affine")
ss1 = "tatggagagaataaaagaactgagagatctaatgtcgcagtcccgcactcgcgagatactcactaagaccactgtggaccatatggccataatcaaaaag".upper()
ss2 = "atggatgtcaatccgactctacttttcctaaaaattccagcgcaaaatgccataagcaccacattcccttatactggagatcctccatacagccatggaa".upper()
print(align(ss1, ss2, score_matrix, gap_open_penalty, gap_extension_penalty))

# Q3
print("report, Q3, global linear")
def pairwise_opt_align(d,score_matrix, gap_penalty):
    n = len(d)
    S = [[0] * (n) for _ in range(n)]

    for i in range(1,n+1): 
        for j in range(1,n+1): 
            alignment1, alignment2, opt_score = needleman_wunsch(d[i], d[j], score_matrix, gap_penalty)
            S[i-1][j-1] = opt_score
    return S

m_linear = pairwise_opt_align(d, score_matrix,gap_penalty)
print(print_matrix(m_linear))

print("report, Q4, global affine")
score_matrix = {
    'A': {'A': 0, 'C': 5, 'G': 2, 'T': 5},
    'C': {'A': 5, 'C': 0, 'G': 5, 'T': 2},
    'G': {'A': 2, 'C': 5, 'G': 0, 'T': 5},
    'T': {'A': 5, 'C': 2, 'G': 5, 'T': 0}
}
gap_open_penalty = 10
gap_extend_penalty = 5

d = {1 : "tatggagagaataaaagaactgagagatctaatgtcgcagtcccgcactcgcgagatactcactaagaccactgtggaccatatggccataatcaaaaag".upper(), 
     2 : "atggatgtcaatccgactctacttttcctaaaaattccagcgcaaaatgccataagcaccacattcccttatactggagatcctccatacagccatggaa".upper(), 
     3 : "tccaaaatggaagactttgtgcgacaatgcttcaatccaatgatcgtcgagcttgcggaaaaggcaatgaaagaatatggggaagatccgaaaatcgaaa".upper(),
     4 : "aaaagcaacaaaaatgaaggcaatactagtagttctgctatatacatttgcaaccgcaaatgcagacacattatgtataggttatcatgcgaacaattca".upper(),
     5 : "atgagtgacatcgaagccatggcgtctcaaggcaccaaacgatcatatgaacaaatggagactggtggggagcgccaggatgccacagaaatcagagcat".upper()}


def pairwise_opt_align_affine(d, score_matrix, gap_open_penalty, gap_extend_penalty):
    n = len(d)
    S = [[0] * (n) for _ in range(n)]

    for i in range(1,n+1): 
        for j in range(1,n+1): 
            alignment, opt_score = align(d[i], d[j], score_matrix, gap_open_penalty, gap_extend_penalty)
            S[i-1][j-1] = opt_score
    return S

m_affine = pairwise_opt_align_affine(d, score_matrix, gap_open_penalty, gap_extend_penalty)
print(print_matrix(m_affine))
