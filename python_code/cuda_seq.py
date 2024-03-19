import blosum as bl
from cuda import cuda, nvrtc

def ASSERT_DRV(err):
    if isinstance(err, cuda.CUresult):
        if err != cuda.CUresult.CUDA_SUCCESS:
            raise RuntimeError("Cuda Error: {}".format(err))
    elif isinstance(err, nvrtc.nvrtcResult):
        if err != nvrtc.nvrtcResult.NVRTC_SUCCESS:
            raise RuntimeError("Nvrtc Error: {}".format(err))
    else:
        raise RuntimeError("Unknown error type: {}".format(err))

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

scoring = bl.BLOSUM(62, default=0)

def translate_codon(codon):
    return CODON_TABLE.get(codon, 'FAIL')

def get_score(scoring, dna, protein, i, j):
    return scoring[translate_codon(dna[i-1:i+2])][protein[j-1]] if translate_codon(dna[i-1:i+2]) != 'FAIL' else -999

def reverse(seq):
    """Returns a reversed string"""
    return seq[::-1]


def complement(seq):
    """Returns a complement DNA sequence"""
    complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    seq_list = list(seq)
    seq_list = [complement_dict[base] for base in seq_list]
    return ''.join(seq_list)


def reverse_complement(seq):
    """"Returns a reverse complement DNA sequence"""
    seq = reverse(seq)
    seq = complement(seq)
    return seq

def matrix_init(dnaSeq, proteinSeq, I, D, C, gep, gop, frameshift_penalty):
    i, j = 0, 0
    N, M = len(dnaSeq), len(proteinSeq)
    TI = [[0 for col in range(M+1)] for row in range(N)]
    TD = [[0 for col in range(M+1)] for row in range(N)]
    TC = [[0 for col in range(M+1)] for row in range(N)]

    for j in range(M+1):
        I[j][0] = float('-inf')
        D[0][j] = D[2][j] = D[3][j] = float('-inf')
        D[1][j] = C[0][j] - gop - gep
        TI[j][0] = float('-inf')
        TD[0][j] = TD[2][j] = TD[3][j] = float('-inf')
        TD[1][j] = 1

    C[0][0] = 0
    TC[0][0] = 0
    for j in range(1, M+1):
        C[0][j] = 0
        C[j][0] = 0
        C[1][j] = max(I[1][j],
                      D[1][j],
                      C[0][j-1] + get_score(scoring, dnaSeq, proteinSeq, 1, j))
        C[2][j] = max(I[2][j],
                      C[0][j-1] + get_score(scoring, dnaSeq, proteinSeq, 2, j) - frameshift_penalty)
        C[3][j] = max(I[3][j],
                      C[1][j-1] + get_score(scoring, dnaSeq, proteinSeq, 3, j) - frameshift_penalty)
        C[4][j] = max(I[4][j],
                      D[4][j],
                      C[1][j-1] + get_score(scoring, dnaSeq, proteinSeq, 4, j),
                      C[2][j-1] + get_score(scoring, dnaSeq, proteinSeq, 4, j) - frameshift_penalty)

        TC[0][j] = 0
        TC[j][0] = 0
        TC[1][j] = -2 if C[1][j] == I[1][j] else -1 if C[1][j] == D[1][j] else 1 if C[1][j] == C[0][j-1] + get_score(scoring, dnaSeq, proteinSeq, 1, j) else 0
        TC[2][j] = -2 if C[2][j] == I[2][j] else 2 if C[2][j] == C[0][j-1] + get_score(scoring, dnaSeq, proteinSeq, 2, j) - frameshift_penalty else 0
        TC[3][j] = -2 if C[3][j] == I[3][j] else 2 if C[3][j] == C[1][j-1] + get_score(scoring, dnaSeq, proteinSeq, 3, j) - frameshift_penalty else 0
        TC[4][j] = -2 if C[4][j] == I[4][j] else -1 if C[4][j] == D[4][j] else 3 if C[4][j] == C[1][j-1] + get_score(scoring, dnaSeq, proteinSeq, 4, j) else 2 if C[4][j] == C[2][j-1] + get_score(scoring, dnaSeq, proteinSeq, 4, j) - frameshift_penalty else 0

    # Turn all negative values in C to 0
    C = [[0 if x < 0 else x for x in row] for row in C]

def modded_three_frame(dnaSeq, proteinSeq, I, D, C, gep, gop, frameshift_penalty):
    for i in range(N):
        for j in range(1, M+1):
            I[i][j] = max(I[i][j-1] - gep, C[i][j-1] - gop - gep)
            if i < 4: continue
            D[i][j] = max(D[i-3][j] - gep, C[i-3][j] - gop - gep)
            C[i][j] = max(I[i][j],
                          D[i][j],
                          C[i-2][j-1] + get_score(scoring, dnaSeq, proteinSeq, i, j) - frameshift_penalty,
                          C[i-3][j-1] + get_score(scoring, dnaSeq, proteinSeq, i, j),
                          C[i-4][j-1] + get_score(scoring, dnaSeq, proteinSeq, i, j) - frameshift_penalty)


if __name__ == '__main__':
    # Initialize CUDA Driver API
    err, = cuda.cuInit(0)
    ASSERT_DRV(err)
    # Retrieve handle for device 0
    err, cuDevice = cuda.cuDeviceGet(0)
    ASSERT_DRV(err)
    # Create context
    err, context = cuda.cuCtxCreate(0, cuDevice)
    ASSERT_DRV(err)

    dnaSeq = 'GGCGTGGCGCAGGCGCAGAGAGGCGCACCGCGCCGGCGCAGGCGCAGAGACACATGCTAGCGCGTCCAGGGGTGGAGGCGTGGCGCAGGCGCAGAGACGCAAGCCTACGGGCGGGGGTTGGGGGGGCGTGTGTTGCAGGAGCAAAGTCGCACGGCGCCGGGCTGGGGCGGGGGGAGGGTGGCGCCGTGCACGCGCAGAAACTCACGTCACGGTGGCGCGGCGCAGAGACGGGTAGAACCTCAGTAATCCGAAAAGCCGGGATCGACCGCCCCTTGCTTGCAGCCGGGCACTACAGGACCC'
    proteinSeq = 'MAGTVLGVGAGVFILALLWVAVLLLCVLLSRASGAARFSVIFLFFGAVIITSVLLLFPRAGEFPAPEVEVKIVDDFFIGRYVLLAFLSAIFLGGLFLVLIHYVLEPIYAKPLHSY'
    N, M = len(dnaSeq), len(proteinSeq)
    gep, gop, frameshift_penalty = 2, 3, 4

    I = [[0 for col in range(M+1)] for row in range(N)]
    D = [[0 for col in range(M+1)] for row in range(N)]
    C = [[0 for col in range(M+1)] for row in range(N)]
    matrix_init(dnaSeq, proteinSeq, I, D, C, gep, gop, frameshift_penalty)
    modded_three_frame(dnaSeq, proteinSeq, I, D, C, gep, gop, frameshift_penalty)
    s1 = max(flatten_list(C))

    I = [[0 for col in range(M+1)] for row in range(N)]
    D = [[0 for col in range(M+1)] for row in range(N)]
    C = [[0 for col in range(M+1)] for row in range(N)]
    dnaSeq = reverse_complement(dnaSeq)
    matrix_init(dnaSeq, proteinSeq, I, D, C, gep, gop, frameshift_penalty)
    modded_three_frame(dnaSeq, proteinSeq, I, D, C, gep, gop, frameshift_penalty)
    s2 = max(flatten_list(C))

    print(f'Score: {max(s1, s2)}')