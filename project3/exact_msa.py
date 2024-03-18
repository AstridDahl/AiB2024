
############################
#        PROJECT 03        #
############################

### Function : sp_exact_3 ###

# Minimize

gap_penalty = 5

score_matrix = {
    'A': {'A': 0, 'C': 5, 'G': 2, 'T': 5},
    'C': {'A': 5, 'C': 0, 'G': 5, 'T': 2},
    'G': {'A': 2, 'C': 5, 'G': 0, 'T': 5},
    'T': {'A': 5, 'C': 2, 'G': 5, 'T': 0}
}

#### code for cost of optimal pairwise (2d) alignment with linear gapcost. 
def empty_matrix2d(r, c):
    m = []
    for i in range(r):
        m.append([])
        for j in range(c):
            m[i].append(None)
    return m

def prepare_matrix2d(r, c, g):
    pm = empty_matrix2d(r, c)
    for i in range(c):
        pm[0][i] = i*g
    for i in range(1, r):
        pm[i][0] = i*g
    return pm

def cost_optimal2k(seq1, seq2, scores, g): # From assignment1.py. 
    r = len(seq1)+1
    c = len(seq2)+1
    fm = prepare_matrix2d(r, c, g)
    for i in range(1, r): # row index
        for j in range(1, c): # col index
            dscore = fm[i-1][j-1] + scores[seq1[i-1]][seq2[j-1]]
            lscore = fm[i][j-1] + g
            uscore = fm[i-1][j] + g
            fm[i][j] = min(dscore, lscore, uscore) # min() since dissimilarity score. 
    return fm[r-1][c-1] # return cost of optimal alignment. r-1 and c-1, since the first row has idx 0 and the first col has idx 0, when matrix is made in python. 

substitution_matrix = {'A': {'A':0, 'T':5, 'G':2, 'C':5},
               'T': {'A':5, 'T':0, 'G':5, 'C':2}, 
               'G': {'A':2, 'T':5, 'G':0, 'C':5}, 
               'C': {'A':5, 'T':2, 'G':5, 'C':0}}

##########################################################

#### Code for cost of optimal seq alignment 3 seq

def empty_matrix3d(r, c, d): 
    m = []
    for i in range(r):
        m.append([]) # make rows.  
        for j in range(c): # make columns. 
            m[i].append([])
            for k in range(d):
                m[i][j].append(None) # make depth.
    return m

def fill_matrix3d(S1, S2, S3, scores, gap_cost): # returns the cost of the optimal alignment of three sequences. 
    r = len(S1) + 1 # n1+1.  
    c = len(S2) + 1 
    d = len(S3) + 1 
    D = empty_matrix3d(r, c, d) # +1 since dimensions of matrix is to be (n1+1)*(n2+1)*(n3+1). 
    # D[r_idx][c_idx][d_idx].
    # Initialization. Correct.
    D[0][0][0] = 0
    # First element in all rows.
    for i in range(r): # i from 0 to n1, since r is not included. The integers returned by range(x) only goes up to x-1.  
        for k in range(d): # k from 0 to n3.  
            D[i][0][k] = cost_optimal2k(S1[0:i], S3[0:k], scores, gap_cost) + (i + k) * gap_cost 
    # First element in all columns. 
    for j in range(c): # j from 0 to n2.
        for k in range(d): 
            D[0][j][k] = cost_optimal2k(S2[0:j], S3[0:k], scores, gap_cost) + (j + k) * gap_cost 
    # All elements in the first depth. 
    for i in range(r): 
        for j in range(c):
            D[i][j][0] = cost_optimal2k(S1[0:i], S2[0:j], scores, gap_cost) + (i + j) * gap_cost # E.g., S1[0:0] gives an empty 
            # string, since the character at the end idx in S1[begin_idx:end_idx] is excluded. I.e., cost of 0 found if empty seqs aligned. 
    # Fill matrix.
    for i in range(len(S1)): # indexing over S1. Ex. len(S1) = 10, then range: 0,1,2...9
        for j in range(len(S2)): # indexing over S2. 
            for k in range(len(S3)): # indexing over S3. 
                cij = scores[S1[i]][S2[j]]
                cik = scores[S1[i]][S3[k]]
                cjk = scores[S2[j]][S3[k]] 

                d1 = D[i][j][k] + cij + cik + cjk # Go directly from 
                d2 = D[i][j][k+1] + cij + 2*gap_cost
                d3 = D[i][j+1][k] + cik + 2*gap_cost
                d4 = D[i+1][j][k] + cjk + 2*gap_cost
                d5 = D[i][j+1][k+1] + 2*gap_cost
                d6 = D[i+1][j][k+1] + 2*gap_cost
                d7 = D[i+1][j+1][k] + 2*gap_cost
                D[i+1][j+1][k+1] = min(d1, d2, d3, d4, d5, d6, d7)
    return D 

def cost_optimal3k(S1, S2, S3, scores, gap_cost): 
    D = fill_matrix3d(S1, S2, S3, scores, gap_cost)
    return D[len(S1)][len(S2)][len(S3)]

############################################################
#### Code for optimal pairwise seq alignment linear gap cost. 
def empty_matrix(r, c): 
    m = []
    for i in range(r):
        m.append([])
        for j in range(c):
            m[i].append(None)
    return m

def prepare_matrix(r, c, gap_cost):
    pm = empty_matrix(r, c)
    for i in range(c):
        pm[0][i] = i*gap_cost
    for i in range(1, r):
        pm[i][0] = i*gap_cost
    return pm

def fill_matrix(S1, S2, scores, gap_cost):
    r = len(S1)+1
    c = len(S2)+1
    fm = prepare_matrix(r, c, gap_cost)
    for i in range(1, r): # row index
        for j in range(1, c): # col index
            dscore = fm[i-1][j-1] + scores[S1[i-1]][S2[j-1]]
            lscore = fm[i][j-1] + gap_cost
            uscore = fm[i-1][j] + gap_cost
            fm[i][j] = min(dscore, lscore, uscore)
    return fm

def get_traceback_arrow(matrix, row, col, match_score, g):
    score_diagonal = matrix[row-1][col-1]
    score_left = matrix[row][col-1]
    score_upper = matrix[row-1][col]

    score_current = matrix[row][col]
    if score_current == score_diagonal + match_score:
        return 'diagonal'
    if score_current == score_left + g:
        return 'left'
    if score_current == score_upper + g:
        return 'up'

def trace_back(seq1, seq2, matrix, score_matrix, g):
    aligned1 = ''
    aligned2 = ''

    row = len(seq1)
    col = len(seq2)

    while row and col > 0:
        base1 = seq1[row-1]
        base2 = seq2[col-1]

        match_score = score_matrix[base1][base2]

        traceback_arrow = get_traceback_arrow(matrix, row, col, match_score, g)

        if traceback_arrow == 'diagonal':
            aligned1 = base1 + aligned1
            aligned2 = base2 + aligned2
            row -= 1
            col -= 1

        elif traceback_arrow == 'up':
            aligned1 = base1 + aligned1
            aligned2 = '-' + aligned2
            row -= 1

        elif traceback_arrow == 'left':
            aligned1 = '-' + aligned1
            aligned2 = base2 + aligned2
            col -= 1 

    while row > 0:
        base1 = seq1[row-1]
        aligned1 = base1 + aligned1
        aligned2 = '-' + aligned2
        row -= 1

    while col > 0:
        base2 = seq2[col-1]
        aligned1 = '-' + aligned1
        aligned2 = base2 + aligned2
        col -= 1  
    return [aligned1, aligned2]

def align(seq1, seq2, scores, g):
    m = fill_matrix(seq1, seq2, scores, g)
    lst = trace_back(seq1, seq2, m, scores, g)
    return lst

###################################################################
#### Code for optimal seq alignment with 3 seq and linear gap cost. 

def MSA3k(S1, S2, S3, scores, gap_cost):
    D = fill_matrix3d(S1, S2, S3, scores, gap_cost)
    aligned1 = ''
    aligned2 = ''
    aligned3 = ''
    i = len(S1) - 1 # index over S1. If e.g., S1 = 'GT', i = 1, 0. 
    j = len(S2) - 1 
    k = len(S3) - 1 
    while i >= 0 and j >= 0 and k >= 0: # > evalueres før 'and'. Indsæt gaps i aligned i while loop body. 
        cij = scores[S1[i]][S2[j]]
        cik = scores[S1[i]][S3[k]]
        cjk = scores[S2[j]][S3[k]]
        # last character of seq1 is seq1[i-1]. last character of S2 is S2[j-1].  
        if D[i+1][j+1][k+1] == D[i][j][k] + cij + cik + cjk:
            aligned1 = S1[i] + aligned1
            aligned2 = S2[j] + aligned2
            aligned3 = S3[k] + aligned3
            i -= 1
            j -= 1
            k -= 1
        elif D[i+1][j+1][k+1] == D[i][j][k+1] + cij + gap_cost: 
            aligned1 = S1[i] + aligned1
            aligned2 = S2[j] + aligned2
            aligned3 = '-' + aligned3
            i -= 1
            j -= 1
        elif D[i+1][j+1][k+1] == D[i][j+1][k] + cik + gap_cost:
            aligned1 = S1[i] + aligned1
            aligned2 = '-' + aligned2
            aligned3 = S3[k] + aligned3
            i -= 1
            k -= 1
        elif D[i+1][j+1][k+1] == D[i+1][j][k] + cjk + gap_cost:
            aligned1 = '-' + aligned1 
            aligned2 = S2[j] + aligned2
            aligned3 = S3[k] + aligned3 
            j -= 1
            k -= 1
        elif D[i+1][j+1][k+1] == D[i][j+1][k+1] + 2 * gap_cost:
            aligned1 = S1[i] + aligned1
            aligned2 = '-' + aligned2
            aligned3 = '-' + aligned3
            i -= 1
        elif D[i+1][j+1][k+1] == D[i+1][j][k+1] + 2 * gap_cost:
            aligned1 = '-' + aligned1 
            aligned2 = S2[j] + aligned2
            aligned3 = '-' + aligned3
            j -= 1
        elif D[i+1][j+1][k+1] == D[i+1][j+1][k] + 2 * gap_cost:
            aligned1 = '-' + aligned1
            aligned2 = '-' + aligned2 
            aligned3 = S3[k] + aligned3
            k -= 1
    if i >= 0 and j >= 0: # k = 0. pairwise seq alignment. linear gap cost. 
        align1, align2 = align(S1[:i+1], S2[:j+1], scores, gap_cost)
        aligned1 = align1 + aligned1
        aligned2 = align2 + aligned2
        aligned3 = '-'*len(align1) + aligned3
        # i and j become <0 here. 
    elif j >= 0 and k >= 0: # i = 0. pairwise seq alignment.
        align2, align3 = align(S2[:i+1], S3[:j+1], scores, gap_cost)
        aligned1 = '-'*len(align2) + aligned1
        aligned2 = align2 + aligned2
        aligned3 = align3 + aligned3
    elif i >= 0 and k >= 0: # j = 0. pairwise seq alignment.
        align1, align3 = align(S1[:i+1], S3[:j+1], scores, gap_cost)
        aligned1 = align1 + aligned1
        aligned2 = '-'*len(align1) + aligned2
        aligned3 = align3 + aligned3
    return [aligned1, aligned2, aligned3]



### Defaults ######
gap_penalty = 5

score_matrix = {
    'A': {'A': 0, 'C': 5, 'G': 2, 'T': 5},
    'C': {'A': 5, 'C': 0, 'G': 5, 'T': 2},
    'G': {'A': 2, 'C': 5, 'G': 0, 'T': 5},
    'T': {'A': 5, 'C': 2, 'G': 5, 'T': 0}
}

##### define a read_fasta function ####
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
    # make the sequences uppercase
    for header, seq in seqs.items():
        seqs[header] = seq.upper()
    return seqs

# ask if the user wants an alignment or just a cost
def ask_user():
    print("Align sequences? (y/n)")
    answer = input()
    if answer == 'y' or answer == 'Y':
        return True
    else:
        return False
#### define a main function ####
def main():
    import sys
    alignment_bool = ask_user()
    if len(sys.argv) != 2:
        print("Usage: python3 exact_msa.py <filename>")
        sys.exit(1)
    file = sys.argv[1]
    seqs = read_fasta(file)
    seqs = list(seqs.values())
    print(cost_optimal3k(seqs[0], seqs[1], seqs[2], score_matrix, gap_penalty))
    if alignment_bool:
        print(MSA3k(seqs[0], seqs[1], seqs[2], score_matrix, gap_penalty))

#### execute the main function ####
if __name__ == '__main__':
    main()