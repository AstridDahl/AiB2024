# 2-approximation algorithm Gusfield. 
import random
substitution_matrix = {'A': {'A':0, 'T':5, 'G':2, 'C':5},
               'T': {'A':5, 'T':0, 'G':5, 'C':2}, 
               'G': {'A':2, 'T':5, 'G':0, 'C':5}, 
               'C': {'A':5, 'T':2, 'G':5, 'C':0}}

# Code to compute cost of optimal pairwise alignment with linear gapcost. 
#### code for cost of optimal pairwise alignment with linear gapcost. 
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

def cost_optimal2k(seq1, seq2, scores, g):
    r = len(seq1)+1
    c = len(seq2)+1
    fm = prepare_matrix2d(r, c, g)
    for i in range(1, r): # row index
        for j in range(1, c): # col index
            dscore = fm[i-1][j-1] + scores[seq1[i-1]][seq2[j-1]]
            lscore = fm[i][j-1] + g
            uscore = fm[i-1][j] + g
            fm[i][j] = min(dscore, lscore, uscore) # min() since dissimilarity score. 
    return fm[r-1][c-1]

# Code to find center string. 
def center_string(Seqs, scores, g): 
    costs = {} 
    for center_string in Seqs: # Try each string as center string at a time.
        Seqs_without_center_string = Seqs[:] 
        Seqs_without_center_string.remove(center_string)
        cost = 0 # sum of all pairwise alignments, where center_string is the first string.  
        for string in Seqs_without_center_string: 
            cost += cost_optimal2k(center_string, string, scores, g)
        costs[center_string] = cost   
    S1 = min(costs, key = costs.get) # min(dictionary, key = dictionary.get) returns the key that has the smallest value. 
    # center_string with the lowest sum of cost. 
    return S1
# print(center_string(['ACGT', 'ATTCT', 'CTCGA'], substitution_matrix, 5)) 

# Code to extend.
def extend(MSA, A): # M is matrix. A is pairwise alignment. 
    length_MSA = len(MSA[0])
    length_A = len(A[0])
    M = []
    for i in range(len(MSA)+1): # Initialize M with one more string than in MSA. 
        M.append('')
    i = 0
    j = 0
    while i < len(MSA[0]) and j < len(A[0]): 
        if MSA[0][i] == '-':
            for r in range(len(MSA)):
                M[r] += MSA[r][i]
            M[r+1] += '-' # extension. 
            i += 1
        elif MSA[0][i] == A[0][j]:
            for r in range(len(MSA)):
                M[r] += MSA[r][i]
            M[r+1] += A[1][j] # extension. 
            i += 1
            j += 1        
        elif A[0][j] == '-':
            for r in range(len(MSA)):
                M[r] += '-'
            M[r+1] += A[1][j]
            j += 1
    # Håndtér, hvis enden af M3 nås før enden af A eller omvendt. 
    while i < len(MSA[0]):
        # Add remaining of M3 to M4 and insert gaps in M4[3]. 
        for r in range(len(MSA)):
            M[r] += MSA[r][i]
        M[r+1] += '-'        
        i += 1
    while j < len(A[0]):
        # Add remaining of A to M4 and insert gaps in strings from M3. 
        for r in range(len(MSA)):
            M[r] += '-'
        M[r+1] += A[1][j]
        j += 1    
    return M

#print(extend(['A--CGT', 'ATTC-T', 'CT-CGA'], ['ACG-T', 'ACGGT'])) 

# Code for optimal pairwise seq alignment linear gap cost. 
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

# Code to make approximate alignment.
def ApproxAlignment(Seqs, scores, g): 
    S1 = center_string(Seqs, scores, g)
    Seqs_without_center_string = Seqs[:]
    Seqs_without_center_string.remove(S1)
    M = [S1] # Initialisering forkert. M skal bestå af de første to sekvenser til start. 
    for i in range(len(Seqs_without_center_string)): 
        A = align(S1, Seqs_without_center_string[i], scores, g)
        M = extend(M, A) 
    return M

# print(ApproxAlignment(['ACGT', 'ATTCT', 'CTCGA', 'ACGGT'], substitution_matrix, 5)) 
# Result should be ['A--CG-T', 'ATTC--T', 'CT-CG-A', 'A--CGGT']

# Code for cost of induced pairwise alignment, linear gap cost. 
def cost_optimal(aligned1, aligned2, scores, g):
    cost = 0
    for i in range(len(aligned1)):
        if aligned1[i] == '-' and aligned2[i] == '-':
            cost = cost
        elif aligned1[i] == '-' or aligned2[i] == '-':
            cost += g
        else: 
            cost += scores[aligned1[i]][aligned2[i]]
    return cost

# Compute cost for MSA linear gap cost. 
def cost_ApproxAlignment(MSA, scores, g):
    cost = 0
    for i in range(len(MSA)):
        for j in range(i+1, len(MSA)): # correct looping. 
            cost += cost_optimal(MSA[i], MSA[j], scores, g)
    return cost

# print(cost_ApproxAlignment(['A--CGT-', 'ATTC-T-', 'CT-CGA-', 'A--CGGT'], substitution_matrix, 5)) 
# print(cost_ApproxAlignment(['A--CG-T', 'ATTC--T', 'CT-CG-A', 'A--CGGT'], substitution_matrix, 5)) 

# Assignment Question 2.
# brca1_bos_taurus = 'ATGGATTTATCTGCGGATCATGTTGAAGAAGTACAAAATGTCCTCAATGCTATGCAGAAAATCTTAGAGTGTCCAATATGTCTGGAGTTGATCAAAGAGCCTGTCTCTACAAAGTGTGACCACATATTTTGCAAATTTTGTATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAATGTCCTTTGTGTAAGAATGA'
# brca1_canis_lupus = 'ATGGATTTATCTGCGGATCGTGTTGAAGAAGTACAAAATGTTCTTAATGCTATGCAGAAAATCTTAGAGTGTCCAATATGTCTGGAGTTGATCAAAGAGCCTGTTTCTACAAAGTGTGATCACATATTTTGCAAATTTTGTATGCTGAAACTTCTCAACCAGAGGAAGGGGCCTTCACAGTGTCCTTTGTGTAAGAACGA'
# brca1_gallus_gallus = 'GCGAAATGTAACACGGTAGAGGTGATCGGGGTGCGTTATACGTGCGTGGTGACCTCGGTCGGTGTTGACGGTGCCTGGGGTTCCTCAGAGTGTTTTGGGGTCTGAAGGATGGACTTGTCAGTGATTGCCATTGGAGACGTGCAAAATGTGCTTTCAGCCATGCAGAAGAACTTGGAGTGTCCAGTCTGTTTAGATGTGAT'
# brca1_homo_sapiens = 'GTACCTTGATTTCGTATTCTGAGAGGCTGCTGCTTAGCGGTAGCCCCTTGGTTTCCGTGGCAACGGAAAAGCGCGGGAATTACAGATAAATTAAAACTGCGACTGCGCGGCGTGAGCTCGCTGAGACTTCCTGGACGGGGGACAGGCTGTGGGGTTTCTCAGATAACTGGGCCCCTGCGCTCAGGAGGCCTTCACCCTCT'
# brca1_macaca_mulatta = 'ATGGATTTATCTGCTGTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCAATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAGATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTGTGTAAGAATGA'

# Center string is brca1_bos_taurus. 
# print(center_string([brca1_bos_taurus, brca1_canis_lupus, brca1_gallus_gallus, brca1_homo_sapiens, brca1_macaca_mulatta], substitution_matrix, 5)) 
# print(brca1_bos_taurus=='ATGGATTTATCTGCGGATCATGTTGAAGAAGTACAAAATGTCCTCAATGCTATGCAGAAAATCTTAGAGTGTCCAATATGTCTGGAGTTGATCAAAGAGCCTGTCTCTACAAAGTGTGACCACATATTTTGCAAATTTTGTATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAATGTCCTTTGTGTAAGAATGA')

# Cost is 3310. 
# print(ApproxAlignment([brca1_bos_taurus, brca1_canis_lupus, brca1_gallus_gallus, brca1_homo_sapiens, brca1_macaca_mulatta], substitution_matrix, 5))
# print(cost_ApproxAlignment(['ATGGATTTATCTGCGGATCATGTTGAAGA-AG-TAC--AA-AAT-G-TCCTCAATGCTATGCA-GAAAATCTTAG--AGTGTCCAAT-ATGTCTGGAGTTGATCAAAGAG-CCT-GTCTCTACAAAGTGTGA-C--CA-C--ATAT-TTTGCAAAT-TTTG--TATGCTGAA-AC-TTCTCAACCA-GAAGAAAGGGCCTTCACAATGTCC--TTTG-TGTAAGAATGA-', 'ATGGATTTATCTGCGGATCGTGTTGAAGA-AG-TAC--AA-AAT-G-TTCTTAATGCTATGCA-GAAAATCTTAG--AGTGTCCAAT-ATGTCTGGAGTTGATCAAAGAG-CCT-GTTTCTACAAAGTGTGA-T--CA-C--ATAT-TTTGCAAAT-TTTG--TATGCTGAA-AC-TTCTCAACCA-GAGGAAGGGGCCTTCACAGTGTCC--TTTG-TGTAAGAACGA-', 'GCGAA---A--TGT-AA-CACGGTAGAGGTGA-T-C--GG-GGT-G--CGTT-ATAC-GTGCGTGGTGACCTCGGTCGGTGT-TGACGGTGCCTGGGGTTCCTCAGAGTGTTTTGGGGTCTGAAGGATG-GA-CT-TGTC--A-GTGATTGCCATT-GGAGA-CGTGCAAAATGTGCTTTCAGCCATGCAG-AAGAA-CTT-GGAGTGTCCAGTCTG-TTTAGATGTGAT', 'GTACCTTGATTT-CGTATTCTG-AGAGGC-TGCTGCTTAGCGGTAGCCCCTTGGT-TTCCGT--GGCAACGGAAA--AGCG-CGGGA-AT-TACAGA-TAAATTAAA-A---CT-GCGACTGCGCGGCGTGAGC-TCG-CTGAGAC-TTCCTGGACGGGGG-ACAGGCTGTG-GG-GTTTC--TCA-GATAACTGGGCCCCTGCGCT-CAG--GAGGCCTTCACCCTCT-', 'ATGGATTTATCTGCTGTTCGCGTTGAAGA-AG-TAC--AA-AAT-G-TCATTAATGCTATGCA-GAAAATCTTAG--AGTGTCCAAT-CTGTCTGGAGTTGATCAAGGAA-CCT-GTCTCCACAAAGTGTGA-C--CA-C--ATAT-TTTGCAGAT-TTTG--CATGCTGAA-AC-TTCTCAACCA-GAAGAAAGGGCCTTCACAGTGTCC--TTTG-TGTAAGAATGA-'], substitution_matrix, 5)) 


### Defaults ######
gap_penalty = 5

score_matrix = {
    'A': {'A': 0, 'C': 5, 'G': 2, 'T': 5},
    'C': {'A': 5, 'C': 0, 'G': 5, 'T': 2},
    'G': {'A': 2, 'C': 5, 'G': 0, 'T': 5},
    'T': {'A': 5, 'C': 2, 'G': 5, 'T': 0}
}

##### define a read_fasta function ####
import random

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
        seqs[header] = ''.join([random.choice(['A', 'C', 'G', 'T']) if char == 'N' else char for char in seq.upper()])
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
    # read the file
    if len(sys.argv) != 2:
        print("Usage: msa_sp_score.py <filename>")
    filename = sys.argv[1]
    seqs = read_fasta(filename)
    seqs = list(seqs.values())
    msa = ApproxAlignment(seqs, score_matrix, gap_penalty)
    # are we aligning or just getting the cost?
    align = ask_user()

    print(cost_ApproxAlignment(msa, score_matrix, gap_penalty))
    if align:
        print(msa)


#### execute the main function ####
if __name__ == '__main__':
    main()