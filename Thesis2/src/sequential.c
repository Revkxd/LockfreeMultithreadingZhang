#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <time.h>

#include "constants.h"
#include "converters.h"
#include "tables.h"
#define PRINTERS 1
#undef PRINTERS

#define MAXROWS 3001
#define MAXCOLS 3001
// #define MAXCOLS 3001 DEFINED IN CONVERTERS
static int gep = 2, gop = 3, frameshift_penalty = 4;
static int ninf = -999;
static int v1, v2, v3, v4, v5;


int letter_to_blosum_index(char letter) {
    switch(letter) {
        case 'A': return 0;
        case 'R': return 1;
        case 'N': return 2;
        case 'D': return 3;
        case 'C': return 4;
        case 'Q': return 5;
        case 'E': return 6;
        case 'G': return 7;
        case 'H': return 8;
        case 'I': return 9;
        case 'L': return 10;
        case 'K': return 11;
        case 'M': return 12;
        case 'F': return 13;
        case 'P': return 14;
        case 'S': return 15;
        case 'T': return 16;
        case 'W': return 17;
        case 'Y': return 18;
        case 'V': return 19;
        case 'B': return 20;
        case 'Z': return 21;
        case 'X': return 22;
        case '*': return 23;
        default: return -1;
    }
}

char* get_codon(char* dnaSequence, int index) {
    char* codon = (char*) malloc(CODON_LENGTH + 1);
    strncpy(codon, dnaSequence + index, CODON_LENGTH);
    codon[CODON_LENGTH] = '\0';
    return codon;
}

char translate_codon(char *codon) {
    int i;
    for(i = 0; i < CODON_TABLE_SIZE; i++) {
        if(strcmp(codon, CODON_TABLE[i].codon) == 0) {
            return CODON_TABLE[i].amino_acid;
        }
    }
    return '0';
}

char get_translated_codon(char* dnaSequence, int index) {
    #if DEBUG
    char* codon = get_codon(dnaSequence, index);
    #else
    char* codon = get_codon(dnaSequence, index-1);
    #endif
    char amino_acid = translate_codon(codon);
    #if DEBUG
    printf("%s, %c\n", codon, amino_acid);
    #endif
    free(codon);
    return amino_acid;
}

int get_score(char amino_acid1, char amino_acid2) {
    if(amino_acid1 == '-' || amino_acid2 == '-')
        return -999;
    int index1 = letter_to_blosum_index(amino_acid1);
    int index2 = letter_to_blosum_index(amino_acid2);
    if(index1 == -1 || index2 == -1)
        return -999;
    return BLOSUM62[index1][index2];
}

void print_matrix(int rows, int cols, int matrix[][MAXCOLS]) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            printf("%4d ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void str_reverse(char* str) {
    int length = strlen(str);

    for (int i = 0, j = length - 1; i < j; ++i, --j) {
        char temp = str[i];
        str[i] = str[j];
        str[j] = temp;
    }
}

void dna_complement(char* dnaSequence) {
    int length = strlen(dnaSequence);

    for (int i = 0; i < length; ++i) {
        switch (dnaSequence[i]) {
            case 'A':
                dnaSequence[i] = 'T';
                break;
            case 'T':
                dnaSequence[i] = 'A';
                break;
            case 'C':
                dnaSequence[i] = 'G';
                break;
            case 'G':
                dnaSequence[i] = 'C';
                break;
        }
    }
}

void reverse_complement(char* dnaSequence) {
    str_reverse(dnaSequence);
    dna_complement(dnaSequence);
    printf("Reverse complement: %s\n", dnaSequence);
}

int max_of_two(int a, int b) {
    return a > b ? a : b;
}

int max_of_three(int a, int b, int c) {
    int temp = a > b ? a : b;
    return temp > c ? temp : c;
}

int max_of_four(int a, int b, int c, int d) {
    int temp1 = a > b ? a : b;
    int temp2 = c > d ? c : d;
    return temp1 > temp2 ? temp1 : temp2;
}

int max_of_five(int a, int b, int c, int d, int e) {
    int temp1 = max_of_three(a, b, c);
    int temp2 = max_of_two(d, e);
    return temp1 > temp2 ? temp1 : temp2;
}

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

    // Print the matrices for debugging
    #ifdef PRINTERS
    printf("I Matrix:\n");
    print_matrix(N, M + 1, I);
    printf("D Matrix:\n");
    print_matrix(N, M + 1, D);
    printf("C Matrix:\n");
    print_matrix(N, M + 1, C);
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

int main(int argc, char *argv[]) {
    // Check if the correct number of command-line arguments is provided
    if (argc != 3) {
        printf("Usage: %s <arg1> <arg2>\n", argv[0]);
        return 1; // Return an error code
    }

    int score;

    double time_taken, start, end;

    char *dna = argv[1];
    char *protein = argv[2];
    printf("DNA Sequence: %s\n", dna);
    printf("Protein Sequence: %s\n", protein);
    start = clock();
    printf("Score: %d\n\n", six_frame(dna, protein));
    end = clock();
    time_taken = (double)(end - start)*1e3 / CLOCKS_PER_SEC;
    printf("Run time taken: %f ms\n\n", time_taken);

    return 0;
}