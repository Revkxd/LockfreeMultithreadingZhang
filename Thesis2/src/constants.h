#pragma once
#define CODON_TABLE_SIZE 64
#define CODON_LENGTH 3
#define NUM_AMINO_ACIDS 26
#define STRING_MAX 1000

typedef struct {
    char codon[4];
    char amino_acid;
} CodonEntry;

typedef char String[STRING_MAX];