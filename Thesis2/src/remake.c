#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <time.h>

#include "maxfunctions.c"
#include "converters.c"

#define PRINTERS 1
// #undef PRINTERS

#define MAXROWS 3001
// #define MAXCOLS 3001 DEFINED IN CONVERTERS
static int gep = 2, gop = 3, frameshift_penalty = 4;
static int ninf = -999;
static int v1, v2, v3, v4, v5;

int I_matrix(int i, int j, int(*I)[MAXCOLS], int(*C)[MAXCOLS]) {
    if (j == 0) {
        return -999;
    }
    else {
        int var1 = I[i][j - 1] - gep;
        int var2 = C[i][j - 1] - gop - gep;
        // printf("i,j: %d %d ||| var1, var2: %d %d ||| max: %d\n",i, j, var1, var2, max_of_two(var1, var2));
        return max_of_two(var1, var2);
    }
}

int D_matrix(int i, int j, int(*D)[MAXCOLS], int(*C)[MAXCOLS]) {
    if (i == 0 || i == 2 || i == 3) {
        return -999;
    }
    else if (i == 1) {
        return C[0][j] - gop - gep;
    }
    else {
        return max_of_two(D[i - 3][j] - gep, C[i - 3][j] - gop - gep);
    }
}

void init_C_matrix(int i, int j, int(*I)[MAXCOLS], int(*D)[MAXCOLS], int(*C)[MAXCOLS], char* dnaSequence, char* proteinSequence, int *maxScore) {
    if ((i >= 0 && j == 0) || (i == 0 && j >= 0)) {
        C[i][j] = 0;
    }
    else if (i == 1 && j > 0) {
        v1 = I[1][j];
        v2 = D[1][j];
        v3 = C[0][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 1));

        C[i][j] = max_of_three(v1, v2, v3);
        if (C[i][j] < 0) {
            C[i][j] = 0;
        }
    }
    else if (i == 2 && j > 0) {
        v1 = I[2][j];
        v2 = C[0][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 2)) - frameshift_penalty; 
        C[i][j] = max_of_two(v1, v2);
        if (C[i][j] < 0) {
            C[i][j] = 0;
        }
    }
    else if (i == 3 && j > 0) {
        v1 = I[3][j];
        v2 = C[1][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 3)) - frameshift_penalty;
        C[i][j] = max_of_two(v1, v2);
        if (C[i][j] < 0) {
            C[i][j] = 0;
        }
    }
}

void fill_C_matrix(int i, int j, int(*I)[MAXCOLS], int(*D)[MAXCOLS], int(*C)[MAXCOLS], char* dnaSequence, char* proteinSequence, int *maxScore) {
    v1 = I[i][j];
    v2 = D[i][j];
    v3 = C[i-2][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, i)) - frameshift_penalty;
    v4 = C[i-3][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, i));
    v5 = C[i-4][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, i)) - frameshift_penalty;
    // printf("Values for %d, %d: %d %d %d %d %d\n", i, j, v1, v2, v3, v4, v5);
    C[i][j] = max_of_five(v1, v2, v3, v4, v5);
}

int modded_three_frame(char* dnaSequence, char* proteinSequence) {
    int N = strlen(dnaSequence), M = strlen(proteinSequence);
    int I[N][MAXCOLS], D[N][MAXCOLS], C[N][MAXCOLS];
    int TI[N][M + 1], TD[N][M + 1], TC[N][M + 1];
    int i,j;
    int maxScore;

    // Anti-Garbage Values
    for(i = 0; i < N; i++) {
        for(j = 0; j < M + 1; j++) {
            I[i][j] = -999;
            D[i][j] = -999;
            C[i][j] = 0;
            TI[i][j] = 0;
            TD[i][j] = 0;
            TC[i][j] = 0;
        }
    }

    // NOTE: INFINITY is a double, just converted to int to save memory
    // Initialization
    C[0][0] = 0;
    TC[0][0] = 0;
    j = 0;
    for(i = 0 ; i < N; i++) {
        I[i][j] = I_matrix(i, j, I, C);
        D[i][j] = D_matrix(i, j, D, C);
        init_C_matrix(i, j, I, D, C, dnaSequence, proteinSequence, &maxScore);
    }

    for(i = 0 ; i < 4; i++) {
        for(j = 1; j < M + 1; j++) {
            I[i][j] = I_matrix(i, j, I, C);
            D[i][j] = D_matrix(i, j, D, C);
            init_C_matrix(i, j, I, D, C, dnaSequence, proteinSequence, &maxScore);
        }
    }

    // Matrix Filling
    for(i = 4; i < N - 1; i++) {
        for(j = 1; j < M + 1; j++) {
            I[i][j] = I_matrix(i, j, I, C);
            D[i][j] = D_matrix(i, j, D, C);
            fill_C_matrix(i, j, I, D, C, dnaSequence, proteinSequence, &maxScore);
        }
    }
    // Traceback Matrix Filling (ignore this first.)
    // for(i = 0; i < N; i++) {
    //     for(j = 0; j < M; j++) {
    //         TI[i][j] = I[i][j] == I[i][j-1] - gep ? 0 : (I[i][j] == C[i][j-1] - gop - gep ? 1 : -999);
    //         TD[i][j] = D[i][j] == D[i-3][j] - gep ? 0 : (D[i][j] == C[i-3][j] - gop - gep ? 1 : -999);
    //         if(C[i][j] == I[i][j])
    //             TC[i][j] = -2;
    //         else if(C[i][j] == D[i][j])
    //             TC[i][j] = -1;
    //         else if(C[i][j] == C[i-2][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, i)) - frameshift_penalty)
    //             TC[i][j] = 2;
    //         else if(C[i][j] == C[i-3][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, i)))
    //             TC[i][j] = 3;
    //         else if(C[i][j] == C[i-4][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, i)) - frameshift_penalty)
    //             TC[i][j] = 4;
    //         else
    //             TC[i][j] = 0;
    //     }
    // }

    // Print the matrices for debugging
    #ifdef PRINTERS
    printf("I Matrix:\n");
    print_matrix(N, M + 1, I);
    printf("D Matrix:\n");
    print_matrix(N, M + 1, D);
    printf("C Matrix:\n");
    print_matrix(N, M + 1, C);
    // printf("TI Matrix:\n");
    // print_matrix(N, M + 1, TI);
    // printf("TD Matrix:\n");
    // print_matrix(N, M + 1, TD);
    // printf("TC Matrix:\n");
    // print_matrix(N, M + 1, TC);
    #endif

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
    int C1[N][M + 1], C2[N][M + 1];
    int max1, max2;
    #ifdef PRINTERS
    printf("First Run:\n");
    #endif
    max1 = modded_three_frame(dnaSequence, proteinSequence);
    
    #ifdef PRINTERS
    printf("Reverse Complement:\n");
    #endif
    reverse_complement(dnaSequence);
    max2 = modded_three_frame(dnaSequence, proteinSequence);

    return max_of_two(max1, max2);
}