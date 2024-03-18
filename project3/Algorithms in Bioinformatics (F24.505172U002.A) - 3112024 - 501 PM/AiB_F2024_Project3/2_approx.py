from Bio import SeqIO

# find center string, S1, such that cost of all pairwise alignments, (S1,S), are minimized
def center_string(seqs, scores, g):
    cost_list = []
    for i in range(len(seqs)):
        cost = 0
        for seq in seqs:
            cost += cost_optimal(seqs[i], seq, scores, g)
        cost_list.append(cost)
    center = cost_list.index(min(cost_list))
    return center

# Storm's code
def extend_msa(M, A):  # M is a multiple alignment represented as a list of columns
    MA = []
    i = 0
    j = 0

    while i < len(M) and j < len(A):

        # Invariant: (1) MA is a valid merge of all columns before column i in M
        # and all columns before column in A, and (2) the first row of M and A up
        # to (but not including) column i and j respectively is the same string
        # if gaps are removed.
        
        if M[i][0] == '-' and A[j][0] == '-':       # Astrid: zero because we are always in the first row, that is the center string
            # Case 1: The next column in MA is column i in M extended with the second symbol
            # in column j in A.
            M[i].append(A[j][1])
            MA.append(M[i])
            i = i + 1
            j = j + 1

        elif M[i][0] == '-' and A[j][0] != '-':
            # Case 2: A[j][0] is a character, so the second symbol in column j in A, A[j][1],
            # must be in the column of MA that is the column in M where the first symbol corresponds
            # to A[j][0]. By the invariant, this column in M is the next column in M, where the first
            # symbol is a character, so we just moved forward in M until we find this column.
            M[i].append('-')
            MA.append(M[i])
            i = i + 1

        elif M[i][0] != '-' and A[j][0] == '-':
            # Case 3: M[i][0] is a character, so column i in M must be in the column of MA that also
            # contains the second symbol from the column in A, where the first symbol is the character
            # corresponding to M[i][0]. By the invariant, this column in A is the next column in A,
            # where the first symbol is a character, so we just add columns from A to MA until we
            # find this column.
            c = ['-']*len(M[i])
            c.append(A[j][1])
            MA.append(c)
            j = j + 1

        elif M[i][0] != '-' and A[j][0] != '-':
            # Case 4: By the invariant the characters M[i][0] and A[j][0] are at the same position
            # in the string spelled by the row of M and A if gaps are removed. The next column in
            # MA is thus column i in M extended with the second symbol in column j in A.
            M[i].append(A[j][1])
            MA.append(M[i])
            i = i + 1
            j = j + 1

    if i < len(M):
        # add the remaining coloumns of M to MA
        while i < len(M):
            MA.append(M[i].append('-'))
            i = i + 1
            
    if j < len(A):
        # add the remaining columns of A to MA
        k = len(MA[-1])
        while j < len(A):
            c = ['-']*(k-1)
            c.append(A[j][1])
            MA.append(c)
            j = j + 1

    return MA

def align_approx(seqs, scores, g):
    center_index = center_string(seqs, scores, g)
    M = seqs[center_index]   
    del seqs[center_index]
    MA = []
    for c in M:
        MA = MA.append([c])        # should be a list of columns (only one row at initialization)
    for seq in seqs:
        seq_aslist = []
        for c in seq:
            seq_aslist = seq_aslist.append([c])
        A =  # optimal pairwise alignment, list of columns (only two rows)
        MA = extend_msa(MA, seq_aslist)   # input is two lists of columns
    
    return MA

#def cost_approx():   # the cost is the sum of pairs

    

# global linear template which we use when finding the center string
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

# Function to read fasta file with one sequence. 
def read_fasta(filename):
    for record in SeqIO.parse(filename, 'fasta'):
        return(record.seq)

seqs_fromslides = ['ACGT','ATTCT','CTCGA']


substitution_matrix = {'A': {'A':0, 'T':5, 'G':2, 'C':5},
               'T': {'A':5, 'T':0, 'G':5, 'C':2}, 
               'G': {'A':2, 'T':5, 'G':0, 'C':5}, 
               'C': {'A':5, 'T':2, 'G':5, 'C':0}}

print(center_string(seqs_fromslides, substitution_matrix, 5))

print(align('AGTGGGTCA', 'AAACTAGCCCG', substitution_matrix, 5))