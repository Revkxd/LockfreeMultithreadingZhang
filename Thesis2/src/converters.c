#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <time.h>

#include "constants.h"
#include "converters.h"
#include "tables.h"

#define DEBUG 1
#undef DEBUG

#define MAXCOLS 3001

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

#ifdef DEBUG
int main() {
    // Test cases for get_score
    printf("Score of A and A: %d\n", get_score('A', 'A'));
    printf("Score of A and R: %d\n", get_score('A', 'R'));
    printf("Score of A and N: %d\n", get_score('A', 'N'));
    printf("Score of A and D: %d\n", get_score('A', 'D'));
    printf("Score of A and C: %d\n", get_score('A', 'C'));

    // Test cases for get_translated_codon
    printf("Translated codon of ATG: %c\n", get_translated_codon("ATGTATG", 0));
    printf("Translated codon of ATG: %c\n", get_translated_codon("ATGTATG", 1));
    printf("Translated codon of ATG: %c\n", get_translated_codon("ATGTATG", 2));
    printf("Translated codon of ATG: %c\n", get_translated_codon("ATGTATG", 3));
    printf("Translated codon of ATG: %c\n", get_translated_codon("ATGTATG", 4));

    return 0;
}
#endif