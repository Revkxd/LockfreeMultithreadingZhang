#pragma once
int letter_to_blosum_index(char letter);

char* get_codon(char* dnaSequence, int index);

char translate_codon(char *codon);

char get_translated_codon(char* dnaSequence, int index);

int get_score(char amino_acid1, char amino_acid2);

void print_matrix(int rows, int cols, int matrix[][cols]);

void str_reverse(char* str);

void dna_complement(char* dnaSequence);

void reverse_complement(char* dnaSequence);
