#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <pthread.h>

#define TABLE_SIZE 500000

typedef struct {
    int i;
    int j;
    int matrix;
    int value;
} entry_t;

// entry_t *hash_table[TABLE_SIZE];

unsigned int hash(int i, int j, int matrix) {
    unsigned int key = i * j * matrix;
    key = (i << 16) | (j & 0xffff) | (matrix << 8);
    return key % TABLE_SIZE;
}

void init_hash_table(entry_t *hash_table[]) {
    for (int i = 0; i < TABLE_SIZE; i++) {
        hash_table[i] = NULL;
    }
}

void print_table(entry_t *hash_table[]) {
    printf("Hash Table\n-------------------\n");
    for (int i = 0; i < TABLE_SIZE; i++) {
        if (hash_table[i] == NULL) {
            printf("\t%i\t---\n", i);
        } else {
            printf("\t%i\t%i\n", i, hash_table[i]->value);
        }
    }
    printf("-------------------\n");
}

void print_table_to_file(FILE *fp, entry_t *hash_table[]) {
    fprintf(fp, "Hash Table\n-------------------\n");
    for (int i = 0; i < TABLE_SIZE; i++) {
        if (hash_table[i] == NULL) {
            fprintf(fp, "\t%i\t---\n", i);
        } else {
            fprintf(fp, "\t%i\t%i\n", i, hash_table[i]->value);
        }
    }
    fprintf(fp, "-------------------\n");
}

void print_entry(entry_t *entry) {
    printf("Entry: %i, %i, %i. Value: %i\n", entry->i, entry->j, entry->matrix, entry->value);
}

void ht_insert(int i, int j, int matrix, int value, entry_t *hash_table[]) {
    // print_entry(entry);
    int index = hash(i, j, matrix);
    entry_t *entry = malloc(sizeof(entry_t));
    entry->i = i;
    entry->j = j;
    entry->matrix = matrix;
    entry->value = value;
    // printf("Inserting at index %i\n", index);
    if (hash_table[index] != NULL) return;
    hash_table[index] = entry;
}

bool ht_search(int i, int j, int matrix, int *value, entry_t *hash_table[]) {
    int index = hash(i, j, matrix);
    if (hash_table[index] != NULL) {
        *value = hash_table[index]->value;
        return true;
    } else {
        return false;
    }
}