import blosum as bl

flatten_list = lambda y:[x for a in y for x in flatten_list(a)] if type(y) is list else [y]
CODON_TABLE = {
    "TTT" : "F", "TTC" : "F", "TTA" : "L", "TTG" : "L",
    "CTT" : "L", "CTC" : "L", "CTA" : "L", "CTG" : "L",
    "ATT" : "I", "ATC" : "I", "ATA" : "I", "ATG" : "M",
    "GTT" : "V", "GTC" : "V", "GTA" : "V", "GTG" : "V",
    "TCT" : "S", "TCC" : "S", "TCA" : "S", "TCG" : "S",
    "CCT" : "P", "CCC" : "P", "CCA" : "P", "CCG" : "P",
    "ACT" : "T", "ACC" : "T", "ACA" : "T", "ACG" : "T",
    "GCT" : "A", "GCC" : "A", "GCA" : "A", "GCG" : "A",
    "TAT" : "Y", "TAC" : "Y", "TAA" : "*", "TAG" : "*",
    "CAT" : "H", "CAC" : "H", "CAA" : "Q", "CAG" : "Q",
    "AAT" : "N", "AAC" : "N", "AAA" : "K", "AAG" : "K",
    "GAT" : "D", "GAC" : "D", "GAA" : "E", "GAG" : "E",
    "TGT" : "C", "TGC" : "C", "TGA" : "*", "TGG" : "W",
    "CGT" : "R", "CGC" : "R", "CGA" : "R", "CGG" : "R",
    "AGT" : "S", "AGC" : "S", "AGA" : "R", "AGG" : "R",
    "GGT" : "G", "GGC" : "G", "GGA" : "G", "GGG" : "G",
}

def translate_codon(codon):
    return CODON_TABLE.get(codon, 'FAIL')

def get_score(scoring, dna, protein, i, j):
    return scoring[translate_codon(dna[i-1:i+2])][protein[j-1]] if translate_codon(dna[i-1:i+2]) != 'FAIL' else -999

def matrix_printer(matrices: list[list]):
    for matrix in matrices:
        for row in matrix:
            print(row)
        print()

def modded_three_frame(dna_input, protein_input):
    N, M = len(dna_input), len(protein_input)
    gep, gop = 2, 3
    frameshift_penalty = 4
    I = [[0 for col in range(M+1)] for row in range(N)]
    D = [[0 for col in range(M+1)] for row in range(N)]
    C = [[0 for col in range(M+1)] for row in range(N)]
    TI = [[0 for col in range(M+1)] for row in range(N)]
    TD = [[0 for col in range(M+1)] for row in range(N)]
    TC = [[0 for col in range(M+1)] for row in range(N)]
    scoring = bl.BLOSUM(62, default=0)

    # Initialization
    for j in range(M+1):
        I[j][0] = float('-inf')
        D[0][j] = D[2][j] = D[3][j] = float('-inf')
        D[1][j] = C[0][j] - gop - gep
        TI[j][0] = float('-inf')
        TD[0][j] = TD[2][j] = TD[3][j] = float('-inf')
        TD[1][j] = 1

    # Note: Placed j-1 for accessing protein_input since index out of bounds error
    C[0][0] = 0
    TC[0][0] = 0
    for j in range(1, M+1):
        C[0][j] = 0
        C[j][0] = 0
        C[1][j] = max(I[1][j],
                      D[1][j],
                      C[0][j-1] + get_score(scoring, dna_input, protein_input, 1, j))
        C[2][j] = max(I[2][j],
                      C[0][j-1] + get_score(scoring, dna_input, protein_input, 2, j) - frameshift_penalty)
        C[3][j] = max(I[3][j],
                      C[1][j-1] + get_score(scoring, dna_input, protein_input, 3, j) - frameshift_penalty)
        C[4][j] = max(I[4][j],
                      D[4][j],
                      C[1][j-1] + get_score(scoring, dna_input, protein_input, 4, j),
                      C[2][j-1] + get_score(scoring, dna_input, protein_input, 4, j) - frameshift_penalty)

        TC[0][j] = 0
        TC[j][0] = 0
        TC[1][j] = -2 if C[1][j] == I[1][j] else -1 if C[1][j] == D[1][j] else 1 if C[1][j] == C[0][j-1] + get_score(scoring, dna_input, protein_input, 1, j) else 0
        TC[2][j] = -2 if C[2][j] == I[2][j] else 2 if C[2][j] == C[0][j-1] + get_score(scoring, dna_input, protein_input, 2, j) - frameshift_penalty else 0
        TC[3][j] = -2 if C[3][j] == I[3][j] else 2 if C[3][j] == C[1][j-1] + get_score(scoring, dna_input, protein_input, 3, j) - frameshift_penalty else 0
        TC[4][j] = -2 if C[4][j] == I[4][j] else -1 if C[4][j] == D[4][j] else 3 if C[4][j] == C[1][j-1] + get_score(scoring, dna_input, protein_input, 4, j) else 2 if C[4][j] == C[2][j-1] + get_score(scoring, dna_input, protein_input, 4, j) - frameshift_penalty else 0

    # Turn all negative values in C to 0
    C = [[0 if x < 0 else x for x in row] for row in C]

    # Matrix filling
    for i in range(N):
        for j in range(1, M+1):
            I[i][j] = max(I[i][j-1] - gep, C[i][j-1] - gop - gep)
            if i < 4: continue
            D[i][j] = max(D[i-3][j] - gep, C[i-3][j] - gop - gep)
            C[i][j] = max(I[i][j],
                          D[i][j],
                          C[i-2][j-1] + get_score(scoring, dna_input, protein_input, i, j) - frameshift_penalty,
                          C[i-3][j-1] + get_score(scoring, dna_input, protein_input, i, j),
                          C[i-4][j-1] + get_score(scoring, dna_input, protein_input, i, j) - frameshift_penalty)

    # Traceback matrix filling
    for i in range(N):
        for j in range(M):
            TI[i][j] = 0 if I[i][j] == I[i][j-1] - gep else 1 if I[i][j] == C[i][j-1] - gop - gep else float('-inf')
            TD[i][j] = 0 if D[i][j] == D[i-3][j] - gep else 1 if D[i][j] == C[i-3][j] - gop - gep else float('-inf')
            if C[i][j] == I[i][j]: TC[i][j] = -2
            elif C[i][j] == D[i][j]: TC[i][j] = -1
            elif C[i][j] == C[i-2][j-1] + get_score(scoring, dna_input, protein_input, i, j) - frameshift_penalty: TC[i][j] = 2
            elif C[i][j] == C[i-3][j-1] + get_score(scoring, dna_input, protein_input, i, j): TC[i][j] = 3
            elif C[i][j] == C[i-4][j-1] + get_score(scoring, dna_input, protein_input, i, j) - frameshift_penalty: TC[i][j] = 4

    matrix_printer([I, D, C, TI, TD, TC])
    return C

def reverse(seq):
    """Returns a reversed string"""
    return seq[::-1]


def complement(seq):
    """Returns a complement DNA sequence"""
    complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    seq_list = list(seq)
    # I can think of 3 ways to do this but first is faster I think ???
    # first list comprehension
    seq_list = [complement_dict[base] for base in seq_list]
    # second complicated lambda
    # seq_list = list(map(lambda base: complement_dict[base], seq_list))
    # third easy for loop
    # for i in range(len(seq_list)):
    #    seq_list[i] = complement_dict[seq_list[i]]
    return ''.join(seq_list)


def reverse_complement(seq):
    """"Returns a reverse complement DNA sequence"""
    seq = reverse(seq)
    seq = complement(seq)
    return seq

def six_frame(dna_input, protein_input):
    C1 = modded_three_frame(dna_input, protein_input)
    print("\n\n")
    C2 = modded_three_frame(reverse_complement(dna_input), protein_input)

    return max(max(flatten_list(C1)), max(flatten_list(C2)))

if __name__ == '__main__':
    dna_inputs = ['ATGCGA', 'ATGCGATACGCTTGA', 'CTTGGTCCGAAT']
    protein_inputs = ['MR', 'MRIR', 'LGPL']
    # thing = bl.BLOSUM(62, default=0)[translate_codon(dna_input[0:3])][protein_input[1]]
    # print(thing)

    for dna_input, protein_input in zip(dna_inputs, protein_inputs):
        ans = six_frame(dna_input, protein_input)
        print(f'Score: {ans}')
        print('-' * 50)