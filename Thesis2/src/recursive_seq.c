#define RECSEQ
#ifdef RECSEQ
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <time.h>

#include "maxfunctions.c"
#include "converters.c"
#include "selfht.c"

#define PRINTERS 1
#undef PRINTERS

int matrixInitialize(char* dnaSequence, char* proteinSequence, int N, int M, int I[][M+1], int D[][M+1], int C[][M+1], int gep, int gop, int frameshift_penalty) {
    int i, j;
    // Anti-Garbage Values
    for(i = 0; i < N; i++) {
        for(j = 0; j < M + 1; j++) {
            I[i][j] = 0;
            D[i][j] = 0;
            C[i][j] = 0;
        }
    }

    // NOTE: INFINITY is a double, just converted to int to save memory
    // Initialization
    for(j = 0; j < M + 1; j++) {
        I[0][j] = -999;
        D[0][j] = -999;
        D[2][j] = -999;
        D[3][j] = -999;
        D[1][j] = C[0][j] - gop - gep;
    }

    C[0][0] = 0;
    for(j = 1; j < M + 1; j++) {
        C[0][j] = 0;
        C[j][0] = 0;
        C[1][j] = max_of_three(I[1][j], D[1][j], C[0][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 1)));
        C[2][j] = max_of_two(I[2][j], C[0][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 2)) - frameshift_penalty);
        C[3][j] = max_of_two(I[3][j], C[1][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 3)) - frameshift_penalty);
        C[4][j] = max_of_four(I[4][j], D[4][j], C[1][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 4)), C[2][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 4)) - frameshift_penalty);
    }

    // Change all negative values in C matrix to 0
    for(i = 0; i < N; i++) {
        for(j = 0; j < M + 1; j++) {
            if(C[i][j] < 0)
                C[i][j] = 0;
        }
    }
}

int calculateI(char* dnaSequence, char* proteinSequence, int i, int j, int gep, int gop, int frameshift_penalty);
int calculateD(char* dnaSequence, char* proteinSequence, int i, int j, int gep, int gop, int frameshift_penalty);
int calculateC(char* dnaSequence, char* proteinSequence, int i, int j, int gep, int gop, int frameshift_penalty);

int calculateI(char* dnaSequence, char* proteinSequence, int i, int j, int gep, int gop, int frameshift_penalty) {
    int score = -999;
    if (ht_search(i, j, 1, &score)) {
        return score;
    }

    // Base case: if j is 0, return 0
    if (i == 0 || j == 0) {
        return 0;
    }

    // Recursive case for I matrix
    score = max_of_two(calculateI(dnaSequence, proteinSequence, i, j - 1, gep, gop, frameshift_penalty) - gep, calculateC(dnaSequence, proteinSequence, i, j - 1, gep, gop, frameshift_penalty) - gop - gep);
    ht_insert(i, j, 1, score);
    return score;
}

int calculateD(char* dnaSequence, char* proteinSequence, int i, int j, int gep, int gop, int frameshift_penalty) {
    int score = -999;
    if(ht_search(i, j, 2, &score)) {
        return score;
    }

    // Base case: if i is less than 4 or j is 0, return 0
    if (i < 4 || j == 0) {
        return 0;
    }

    // Recursive case for D matrix
    score = max_of_two(calculateD(dnaSequence, proteinSequence, i - 3, j, gep, gop, frameshift_penalty) - gep, calculateC(dnaSequence, proteinSequence, i - 3, j, gep, gop, frameshift_penalty) - gop - gep);
    ht_insert(i, j, 2, score);
    return score;
}

int calculateC(char* dnaSequence, char* proteinSequence, int i, int j, int gep, int gop, int frameshift_penalty) {
    int score = -999;
    if (ht_search(i, j, 3, &score)) {
        return score;
    }

    // Base case: if i or j is 0, return 0
    if (i == 0 || j == 0) {
        return 0;
    }

    // Recursive cases for C matrix
    int I_val = calculateI(dnaSequence, proteinSequence, i, j, gep, gop, frameshift_penalty);
    int D_val = calculateD(dnaSequence, proteinSequence, i, j, gep, gop, frameshift_penalty);

    // Calculate C[i][j] based on I[i][j], D[i][j], and C[i][j]
    // Adjust this part based on your actual calculation logic
    score = max_of_five(
        I_val,
        D_val,
        calculateC(dnaSequence, proteinSequence, i - 2, j - 1, gep, gop, frameshift_penalty) + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, i)) - frameshift_penalty,
        calculateC(dnaSequence, proteinSequence, i - 3, j - 1, gep, gop, frameshift_penalty) + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, i)),
        calculateC(dnaSequence, proteinSequence, i - 4, j - 1, gep, gop, frameshift_penalty) + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, i)) - frameshift_penalty
    );
    ht_insert(i, j, 3, score);
    return score;
}

int modded_three_frame(char* dnaSequence, char* proteinSequence, int N, int M, int I[][M+1], int D[][M+1], int C[][M+1], int gep, int gop, int frameshift_penalty) {
    int i, j;

    matrixInitialize(dnaSequence, proteinSequence, N, M, I, D, C, gep, gop, frameshift_penalty);

    // Matrix Filling
    for(i = 0; i < N; i++) {
        for(j = 1; j < M + 1; j++) {
            I[i][j] = calculateI(dnaSequence, proteinSequence, i, j, gep, gop, frameshift_penalty);
            // I[i][j] = max_of_two(I[i][j-1] - gep, C[i][j-1] - gop - gep);
            if (i < 4)
                continue;
            D[i][j] = calculateD(dnaSequence, proteinSequence, i, j, gep, gop, frameshift_penalty);
            C[i][j] = calculateC(dnaSequence, proteinSequence, i, j, gep, gop, frameshift_penalty);
        }
    }

    int max_val = -999;
    for(i = 0; i < N; i++) {
        for(j = 0; j < M + 1; j++) {
            max_val = C[i][j] > max_val ? C[i][j] : max_val;
        }
    }
    
    return max_val;
}

int six_frame(char* dnaSequence, char* proteinSequence) {
    int N = strlen(dnaSequence), M = strlen(proteinSequence);
    int gep = 2, gop = 3, frameshift_penalty = 4;
    int I[N][M + 1], D[N][M + 1], C[N][M+1];
    int TI[N][M + 1], TD[N][M + 1], TC[N][M + 1];
    int max1, max2;
    #ifdef PRINTERS
    printf("First Run:\n");
    #endif
    max1 = modded_three_frame(dnaSequence, proteinSequence, N, M, I, D, C, gep, gop, frameshift_penalty);

    #ifdef PRINTERS
    printf("Reverse Complement:\n");
    #endif
    reverse_complement(dnaSequence);
    max2 = modded_three_frame(dnaSequence, proteinSequence, N, M, I, D, C, gep, gop, frameshift_penalty);

    return max_of_two(max1, max2);
}

int main() {
    String dnaSequences[] = {"ATTGACAACCGCGTCCGCCGC","ATTGACAACCGCGTCCGCCGCCGCTTCAAGGGCCAGTACTTGATGCCCAACATTGGCTACGGCTCCAACAAGCGCACCCGCCACATGTTGCCCACCGGCT", "GCTACGTCCGCTCCTCCATGTCCTTGTCCGGCTACATGCCCCCCTTGTGCGACCCCAAGGACGGCCACTTGTTGTTGGACGGCGGCTACGTCAACAACT", "GAGCCCACCTCCGAGATTTTGCAGAACCCCGCCCGCGTCTTGCGCCAGCAGTTGAAGGTCTTGTCCGTCATTGACGGCCAGTCCTACGAGCCCTTGAAGG", "CCCGGCGCCGGCTCCGGCCACGGCCACGGCCCCAACGGCGGCTCCAACTCCTCCTCCTGCACCCCCCCCTCCTCCAACCCCCACATTACCGGCTACGTCG"};
    String proteinSequences[] = {"IDNRVR","IDNRVRRRFKGQYLMPNIGYGSNKRTRHMLPTGF", "RYVRSSMSLSGYMPPLCDPKDGHLLLDGGYVNNL", "EPTSEILQNPARVLRQQLKVLSVIDGQSYEPLKD", "PGAGSGHGHGPNGGSNSSSCTPPSSNPHITGYVD"};
    int i;
    double time_taken, start, end;
    for(i = 0; i < 5; i++) {
        init_hash_table();
        printf("DNA Sequence: %s\n", dnaSequences[i]);
        printf("Protein Sequence: %s\n", proteinSequences[i]);
        start = clock();
        printf("Score: %d\n\n", six_frame(dnaSequences[i], proteinSequences[i]));
        end = clock();
        time_taken = (double)(end - start)*1e3 / CLOCKS_PER_SEC;
        printf("Run %d time taken: %f ms\n\n", i, time_taken);
        // break;
    }

    // String dnaSeq = "GGCGTGGCGCAGGCGCAGAGAGGCGCACCGCGCCGGCGCAGGCGCAGAGACACATGCTAGCGCGTCCAGGGGTGGAGGCGTGGCGCAGGCGCAGAGACGCAAGCCTACGGGCGGGGGTTGGGGGGGCGTGTGTTGCAGGAGCAAAGTCGCACGGCGCCGGGCTGGGGCGGGGGGAGGGTGGCGCCGTGCACGCGCAGAAACTCACGTCACGGTGGCGCGGCGCAGAGACGGGTAGAACCTCAGTAATCCGAAAAGCCGGGATCGACCGCCCCTTGCTTGCAGCCGGGCACTACAGGACCC";
    // String protSeq = "PROHISARGVALARGVALSERPROARGGLYALAALAALASERALASERLEUCYSTHRILEALAGLNVALPROTHRSERALAPRORGLYVALARGMETPROALAPRONPROALAHISASNVALLEUVALSERALACYSARGGLYPROTHRPROPROPROSERHISARGGLYTHRCYSALASERLEUSERALAVAPRORARGARGVALSERALAHILEUGLYVALILEARGLEUPHEGLYPROSERTRPARGGLYTHRASNVALGLYPROCYSPROGLY";
    
    // char dnaSeq[STRING_MAX];
    // strcpy(dnaSeq,
    //       "ATTGACAACCGCGTCCGCCGC"
    //       );
    // char protSeq[STRING_MAX];
    // strcpy(protSeq,
    //       "IDNRVR"
    //       );
    // start = clock();
    // printf("Score: %d\n\n", six_frame(dnaSeq, protSeq));
    // end = clock();
    // time_taken = (double)(end - start)*1e3 / CLOCKS_PER_SEC;
    // printf("Run %d time taken: %f ms\n\n", i, time_taken);
    return 0;
}
#endif