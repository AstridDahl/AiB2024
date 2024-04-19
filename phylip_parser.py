
def parse_phylip(file):
    '''Takes phylip file with distance matrix and returns dictionary with distance matrix.'''
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

print(parse_phylip('example_slide4.phy'))