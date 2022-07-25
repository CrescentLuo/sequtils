from os import scandir
import numpy as np

'''
A={a_1,a_2,...a_n}
B={b_1,b_2,...b_m}

H_ij = Max{
    H_i-1,j-1 + s(a_i, b_j)
    max k>=1 {H_i-k,j - W_k} (1<=i<=n, 1<=j<=m)
    max l>=1 {H_i,j-l - W_l} (i<=i<=n, 1<=j<=m)
    0
}
'''
def scoring_matrix(a,b):
    return None
    #return scoring_matrix[iupac_dict[a], iupac_dict[b]]   

def smith_waterman(seq1, seq2):
    iupac_dict = {
        'T': 0,
        'C': 1,
        'A': 2,
        'G': 3,
        '-': 4,
        'R': 5,
        'Y': 6,
        'W': 7,
        'S': 8,
        'M': 9,
        'K': 10,
        'H': 11,
        'B': 12,
        'V': 13,
        'D': 14,
        'N': 15,
        '?': 16}

    match_score = np.array([
        [36, 0, 0, 0, 9, 0, 18, 18, 0, 0, 18, 12, 12, 0, 12, 9, 0],
        [0, 36, 0, 0, 9, 0, 18, 0, 18, 18, 0, 12, 12, 12, 0, 9, 0],
        [0, 0, 36, 0, 9, 18, 0, 18, 0, 18, 0, 12, 0, 12, 12, 9, 0],
        [0, 0, 0, 36, 9, 18, 0, 0, 18, 0, 18, 0, 12, 12, 12, 9, 0],
        [9, 9, 9, 9, 36, 18, 18, 18, 18, 18, 18, 27, 27, 27, 27, 36, 0],
        [0, 0, 18, 18, 18, 36, 0, 9, 9, 9, 9, 6, 6, 12, 12, 18, 0],
        [18, 18, 0, 0, 18, 0, 36, 9, 9, 9, 9, 12, 12, 6, 6, 18, 0],
        [18, 0, 18, 0, 18, 9, 9, 36, 0, 9, 9, 12, 6, 6, 12, 18, 0],
        [0, 18, 0, 18, 18, 9, 9, 0, 36, 9, 9, 6, 12, 12, 6, 18, 0],
        [0, 18, 18, 0, 18, 9, 9, 9, 9, 36, 0, 12, 6, 12, 6, 18, 0],
        [18, 0, 0, 18, 18, 9, 9, 9, 9, 0, 36, 6, 12, 6, 12, 18, 0],
        [12, 12, 12, 0, 27, 6, 12, 12, 6, 12, 6, 36, 8, 8, 8, 27, 0],
        [12, 12, 0, 12, 27, 6, 12, 6, 12, 6, 12, 8, 36, 8, 8, 27, 0],
        [0, 12, 12, 12, 27, 12, 6, 6, 12, 12, 6, 8, 8, 36, 8, 27, 0],
        [12, 0, 12, 12, 27, 12, 6, 12, 6, 6, 12, 8, 8, 8, 36, 27, 0],
        [9, 9, 9, 9, 36, 18, 18, 18, 18, 18, 18, 27, 27, 27, 27, 36, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]], dtype=object)

    seq1_len = len(seq1)
    seq2_len = len(seq2)
    score_matrix = np.zeros([seq1_len+1, seq2_len+1])
    for i, a in enumerate(seq1,1):
        for j, b in enumerate(seq2,1):
            score_match = score_matrix[i-1,j-1] + match_score(a,b)
            score_gap_a = score_matrix[i-1,j] + match_score(a, '-')
            score_gap_b = score_matrix[i, j-1] + match_score('-', b)
            score_matrix = max(score_match, score_gap_a, score_gap_b, 0)

    x, y = np.unravel_index(score_matrix.argmax(), score_matrix.shape)
    path_a = ''
    path_b = ''
    pos_score = score_matrix[x,y]
    path_a += seq1[x-1]
    path_b += seq2[y-1]

    while pos_score > 0:

        pos_score = score_matrix[x,y]
        if (score_matrix[x-1, y-1] >= score_matrix[x,y-1]) and (score_matrix[x-1, y-1] >= score_matrix[x-1, y]):
            x -= 1
            y -= 1
            path_a += seq1[x-1]
            path_b += seq2[y-1]

        if (score_matrix[x-1, y] >= score_matrix[x-1, y-1]) and (score_matrix[x-1, y] >= score_matrix[x, y-1]):
            x -= 1
            path_a += seq1[x-1]
            path_b += '-'

        if (score_matrix[x, y-1] >= score_matrix[x-1, y-1]) and (score_matrix[x, y-1] >= score_matrix[x-1, y]):
            y -= 1
            path_a += '-'
            path_b += seq2[y-1]

    print(path_a[::-1])
    print(path_b[::-1])    

def compare_matrix(ma_x, ma_y):
    n, m = ma_x.shape[0], ma_y.shape[0]
    pcc_matrix = np.zeros((n,m))
    for i in range(n):
        for j in range(m):
            if j > i:
                continue
            pcc_matrix[i,j] = np.corrcoef(ma_x[i,:], ma_y[j,:])[0,1]
    return pcc_matrix

