# affine gap cost. g(k) = a + b(k-1), where a is opening cost and b is extension cost. 
# Here, g(k) = 5+5k, i.e., the opening cost a = 10 and the extension cost b = 5.
# dissimilarity score. 

import numpy as np

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
        pm[0][j] = a+b*(j-1) # In book, -a and -b, since a similarity 
        # score is to be maximized. Here +a and +b since a dissimilarity 
        # score is to be minimized.
    for i in range(1, r):
        pm[i][0] = a+b*(i-1)
    return pm
# print(prepare_matrix_D(4,4, 10, 5))

def prepare_matrix_E(r, c):
    pm = empty_matrix(r, c)
    for i in range(c):
        pm[0][i] = np.inf 
    return pm
# print(prepare_matrix_E(4,4)) 

def prepare_matrix_F(r, c):
    pm = empty_matrix(r, c)
    for i in range(r):
        pm[i][0] = np.inf
    return pm 
# print(prepare_matrix_F(4,4))

substitution_matrix = {'A': {'A':0, 'T':5, 'G':2, 'C':5},
               'T': {'A':5, 'T':0, 'G':5, 'C':2}, 
               'G': {'A':2, 'T':5, 'G':0, 'C':5}, 
               'C': {'A':5, 'T':2, 'G':5, 'C':0}}

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

def fill_matrices(seq1, seq2, scores, a, b): # Giver korrekte matricer. Tjekket med Johans eksempel. 
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
    # Behøver egentlig kun matrix D for at lave backtrack...
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
    matrices = fill_matrices(seq1, seq2, scores, a, b)
    alignment = traceback(seq1, seq2, matrices, scores, a, b)
    return alignment

# Example
seq1 = 'aggt'
seq2 = 'acta'
# print(fill_matrices(seq1, seq2, substitution_matrix, 10, 5))

# Case 1: score of 24. alignment works. 
seq1 = 'acgtgtcaacgt' 
seq2 = 'acgtcgtagcta'
# print(compute_cost(seq1, seq2, substitution_matrix, 10, 5)) 
# print(align(seq1, seq2, substitution_matrix, 10, 5))

# Case 2: score of 22. works. 
seq1 = 'aataat' # 'aataat'
seq2 = 'aagg' # 'aagg'
# print(compute_cost(seq1, seq2, substitution_matrix, 10, 5))
# matrices = fill_matrices(seq1, seq2, substitution_matrix, 10, 5)
# print(align(seq1, seq2, substitution_matrix, 10, 5))

# Case 3: score of 29. Works. 
seq1 = 'tccagaga'
seq2 = 'tcgat'
# print(compute_cost(seq1, seq2, substitution_matrix, 10, 5))
print(fill_matrices(seq1, seq2, substitution_matrix, 10, 5))
print(align(seq1, seq2, substitution_matrix, 10, 5))

# Case 4: score of 395. Works. 
seq1 = 'ggcctaaaggcgccggtctttcgtaccccaaaatctcggcattttaagataagtgagtgttgcgttacactagcgatctaccgcgtcttatacttaagcgtatgcccagatctgactaatcgtgcccccggattagacgggcttgatgggaaagaacagctcgtctgtttacgtataaacagaatcgcctgggttcgc'
seq2 = 'gggctaaaggttagggtctttcacactaaagagtggtgcgtatcgtggctaatgtaccgcttctggtatcgtggcttacggccagacctacaagtactagacctgagaactaatcttgtcgagccttccattgagggtaatgggagagaacatcgagtcagaagttattcttgtttacgtagaatcgcctgggtccgc'
# print(compute_cost(seq1, seq2, substitution_matrix, 10, 5))
# print(align(seq1, seq2, substitution_matrix, 10, 5))
# print('ggcctaaaggcgccggtctttcgtaccccaaaatctcggcattttaagataagtgagtgttgcgttacactagcgatctaccgcgtcttatacttaagcgtatgcccagatctgactaatcgtgcccccggattagacgggcttgatgggaaagaacagctcgtc------tgtttacgtataaacagaatcgcctgggttcgc' == 'ggcctaaaggcgccggtctttcgtaccccaaaatctcggcattttaagataagtgagtgttgcgttacactagcgatctaccgcgtcttatacttaagcgtatgcccagatctgactaatcgtgcccccggattagacgggcttgatgggaaagaacagctcgtc------tgtttacgtataaacagaatcgcctgggttcgc')
# print('gggctaaaggttagggtctttcacactaaagagtggt-gcgtatcgtggctaatgtaccgcttctggtatc-gtggcttacggc--cagacctacaagtactagacctga--gaactaatcttgtcgagccttccattgagggtaatgggagagaacatcgagtcagaagttattcttgtttacgtagaatcgcctgggtccgc' == 'gggctaaaggttagggtctttcacactaaagagtggt-gcgtatcgtggctaatgtaccgcttctggtatc-gtggcttacggc--cagacctacaagtactagacctga--gaactaatcttgtcgagccttccattgagggtaatgggagagaacatcgagtcagaagttattcttgtttacgtagaatcgcctgggtccgc')
