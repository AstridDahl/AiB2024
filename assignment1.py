from Bio import SeqIO

def empty_matrix(r, c):
    m = []
    for i in range(r):
        m.append([])
        for j in range(c):
            m[i].append(None)
    return m

def prepare_matrix(r, c, g):
    pm = empty_matrix(r, c)
    for i in range(c):
        pm[0][i] = i*g
    for i in range(1, r):
        pm[i][0] = i*g
    return pm

def fill_matrix(seq1, seq2, scores, g):
    r = len(seq1)+1
    c = len(seq2)+1
    fm = prepare_matrix(r, c, g)
    for i in range(1, r): # row index
        for j in range(1, c): # col index
            dscore = fm[i-1][j-1] + scores[seq1[i-1]][seq2[j-1]]
            lscore = fm[i][j-1] + g
            uscore = fm[i-1][j] + g
            fm[i][j] = max(dscore, lscore, uscore)
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

def cost_optimal(seq1, seq2, scores, g):
    align1, align2 = align(seq1, seq2, scores, g)
    cost = 0
    for i in range(len(align1)):
        if align1[i] == '-' or align2[i] == '-':
            cost += g
        else: 
            cost += scores[align1[i]][align2[i]]
    return cost

substitution_matrix = {'A': {'A':10, 'T':2, 'G':5, 'C':2},
               'T': {'A':2, 'T':10, 'G':2, 'C':5}, 
               'G': {'A':5, 'T':2, 'G':10, 'C':2}, 
               'C': {'A':2, 'T':5, 'G':2, 'C':10}}

# Function to read fasta file with one sequence. 
def read_fasta(filename):
    for record in SeqIO.parse(filename, 'fasta'):
        return(record.seq)
    
seq1 = read_fasta('seq1.fasta')
seq2 = read_fasta('seq2.fasta')

# Question 2
print(cost_optimal(seq1, seq2, substitution_matrix, -5))

# Question 3
alignment = align(seq1, seq2, substitution_matrix, -5)
for s in alignment:
    print(s)

# Question 4

# Tests. 
A = 'TCCAGAGA'
B = 'TCGAT'
print(cost_optimal(A, B, substitution_matrix, -5))
alignment = align(A, B, substitution_matrix, -5)
for s in alignment:
    print(s)

C = 'CGTGTCAAGTCT'
D = 'ACGTCGTAGCTAGG'
print(cost_optimal(C, D, substitution_matrix, -5))
alignment = align(C, D, substitution_matrix, -5)
for s in alignment:
    print(s)
