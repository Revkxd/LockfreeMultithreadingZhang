def edit_distance(seq1, seq2):
    size = len(seq1)
    scoring = [[0 for _ in range(size + 1)] for _ in range(size + 1)]
    for i in range(size + 1):
        scoring[i][0] = i
        scoring[0][i] = i
    
    for i in range(1, size + 1):
        for j in range(1, size + 1):
            if seq1[i - 1] == seq2[j - 1]:
                scoring[i][j] = scoring[i - 1][j - 1]
            else:
                scoring[i][j] = min(scoring[i - 1][j - 1], scoring[i - 1][j], scoring[i][j - 1]) + 1

    return scoring, scoring[size][size]

if __name__ == '__main__':
    seq1 = 'GATTACA'
    seq2 = 'GCATGCU'

    scoring, distance = edit_distance(seq1, seq2)
    for row in scoring:
        print(row)
    print(f'Edit distance: {distance}')