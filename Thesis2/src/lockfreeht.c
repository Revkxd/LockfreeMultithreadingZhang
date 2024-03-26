#define RECURSIVE
#ifdef RECURSIVE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <atomicdefs.h>

#define TABLE_SIZE 500000

typedef struct {
    int i;
    int j;
    int matrix;
    int value;
} entry_t;

typedef struct ht_ent {
    struct ht_ent *next;
    int i;
    int j;
    int matrix;
    int value;
} ht_entry;

static ht_entry *hash_table[TABLE_SIZE];

unsigned int hash(int i, int j, int matrix) {
    unsigned int key = i * j * matrix;
    key = (i << 16) | (j & 0xffff) | (matrix << 8);
    return key % TABLE_SIZE;
}

void init_hash_table() {
    for (int i = 0; i < TABLE_SIZE; i++) {
        hash_table[i] = NULL;
    }
}

void print_table() {
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

void print_table_to_file(FILE *fp) {
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

void print_entry(ht_entry *entry) {
    printf("Entry: %i, %i, %i\n", entry->i, entry->j, entry->matrix);
}

void ht_insert(int i, int j, int matrix, int value) {
    // print_entry(entry);
    unsigned int index = hash(i, j, matrix);
    ht_entry *entry = (ht_entry *)malloc(sizeof(ht_entry));
    entry->value = value;
    // printf("Inserting at index %i\n", index);
    ht_entry *ent, *oldent, *newent = NULL;
    ht_entry *inserted_entry = NULL;
    if (hash_table[index] != NULL) return;
    hash_table[index] = entry;

    do {
        while (ent && (i == ent->i && j == ent->j && matrix == ent->matrix)) {
            oldent = ent;
            ent = ent->next;
        }
    } while (hash_table[index] != NULL);
}

bool ht_search(int i, int j, int matrix, int *value) {
    int index = hash(i, j, matrix);
    if (hash_table[index] != NULL) {
        *value = hash_table[index]->value;
        return true;
    } else {
        return false;
    }
}

#endif