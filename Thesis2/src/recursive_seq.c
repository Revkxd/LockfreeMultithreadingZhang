#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <time.h>

#include "maxfunctions.c"
#include "converters.c"
#include "splitlfht.c"

int calculateI(char* dnaSequence, char* proteinSequence, int i, int j, int gep, int gop, int frameshift_penalty);
int calculateD(char* dnaSequence, char* proteinSequence, int i, int j, int gep, int gop, int frameshift_penalty);
int calculateC(char* dnaSequence, char* proteinSequence, int i, int j, int gep, int gop, int frameshift_penalty);

int matrixInitialize(char* dna, char* protein, int N, int M, int I[][M+1], int D[][M+1], int C[][M+1], int TI[][M+1], int TD[][M+1], int TC[][M+1], int gep, int gop, int frameshift_penalty) {
    int j;
    for(j = 0; j < M + 1; j++) {
        ht_insert(0, j, 1, -999);
        ht_insert(0, j, 2, -999);
        ht_insert(2, j, 2, -999);
        ht_insert(3, j, 2, -999);
        ht_insert(1, j, 2, calculateC(dna, protein, 0, j, gep, gop, frameshift_penalty) - gop - gep);
    }

    for(j = 1; j < M + 1; j++) {
        ht_insert(0, j, 3, 0);
        ht_insert(j, 0, 3, 0);
        ht_insert(1, j, 3, max_of_three(calculateI(dna, protein, 1, j, gep, gop, frameshift_penalty),
                                        calculateD(dna, protein, 1, j, gep, gop, frameshift_penalty),
                                        calculateC(dna, protein, 0, j - 1, gep, gop, frameshift_penalty) + get_score(protein[j - 1], get_translated_codon(dna, 1))));
        ht_insert(2, j, 3, max_of_two(calculateI(dna, protein, 2, j, gep, gop, frameshift_penalty),
                                    calculateC(dna, protein, 0, j - 1, gep, gop, frameshift_penalty) + get_score(protein[j - 1], get_translated_codon(dna, 2)) - frameshift_penalty));
        ht_insert(3, j, 3, max_of_two(calculateI(dna, protein, 3, j, gep, gop, frameshift_penalty),
                                    calculateC(dna, protein, 1, j - 1, gep, gop, frameshift_penalty) + get_score(protein[j - 1], get_translated_codon(dna, 3)) - frameshift_penalty));
        ht_insert(4, j, 3, max_of_four(calculateI(dna, protein, 4, j, gep, gop, frameshift_penalty),
                                    calculateD(dna, protein, 4, j, gep, gop, frameshift_penalty),
                                    calculateC(dna, protein, 1, j - 1, gep, gop, frameshift_penalty) + get_score(protein[j - 1], get_translated_codon(dna, 4)),
                                    calculateC(dna, protein, 2, j - 1, gep, gop, frameshift_penalty) + get_score(protein[j - 1], get_translated_codon(dna, 4)) - frameshift_penalty));
    }
}

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

int modded_three_frame(char* dnaSequence, char* proteinSequence, int N, int M, int I[][M+1], int D[][M+1], int C[][M+1], int TI[][M+1], int TD[][M+1], int TC[][M+1], int gep, int gop, int frameshift_penalty) {
    int i, j;

    matrixInitialize(dnaSequence, proteinSequence, N, M, I, D, C, TI, TD, TC, gep, gop, frameshift_penalty);

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

    // Traceback Matrix Filling
    for(i = 0; i < N; i++) {
        for(j = 0; j < M; j++) {
            TI[i][j] = I[i][j] == I[i][j-1] - gep ? 0 : (I[i][j] == C[i][j-1] - gop - gep ? 1 : -999);
            TD[i][j] = D[i][j] == D[i-3][j] - gep ? 0 : (D[i][j] == C[i-3][j] - gop - gep ? 1 : -999);
            if(C[i][j] == I[i][j])
                TC[i][j] = -2;
            else if(C[i][j] == D[i][j])
                TC[i][j] = -1;
            else if(C[i][j] == C[i-2][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, i)) - frameshift_penalty)
                TC[i][j] = 2;
            else if(C[i][j] == C[i-3][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, i)))
                TC[i][j] = 3;
            else if(C[i][j] == C[i-4][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, i)) - frameshift_penalty)
                TC[i][j] = 4;
            else
                TC[i][j] = 0;
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
    max1 = modded_three_frame(dnaSequence, proteinSequence, N, M, I, D, C, TI, TD, TC, gep, gop, frameshift_penalty);

    #ifdef PRINTERS
    printf("Reverse Complement:\n");
    #endif
    reverse_complement(dnaSequence);
    max2 = modded_three_frame(dnaSequence, proteinSequence, N, M, I, D, C, TI, TD, TC, gep, gop, frameshift_penalty);

    return max_of_two(max1, max2);
}