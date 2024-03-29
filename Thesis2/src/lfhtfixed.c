#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdatomic.h>
#include "atomicdefs.h"

#define TABLE_SIZE 200000

typedef struct ht_ent {
    struct ht_ent *next;
    int i;
    int j;
    int matrix;
    int score;
} ht_entry;

static ht_entry *hash_table[TABLE_SIZE];

void init_hash_table() {
    for (int i = 0; i < TABLE_SIZE; i++) {
        hash_table[i] = NULL;
    }
}

unsigned int hash(int i, int j, int matrix) {
    unsigned int key = (i << 16) | (j & 0xffff) | (matrix << 8);
    return key % TABLE_SIZE;
}

int keymatch(ht_entry *entry, int i, int j, int matrix) {
    return entry->i == i && entry->j == j && entry->matrix == matrix;
}

void ht_insert(int i, int j, int matrix, int score) {
    unsigned int key = hash(i, j, matrix);
    ht_entry *newent = (ht_entry *)malloc(sizeof(ht_entry));
    if (!newent) {
        fprintf(stderr, "Memory allocation failed\n");
        return;
    }

    newent->i = i;
    newent->j = j;
    newent->matrix = matrix;
    newent->score = score;
    newent->next = NULL;

    ht_entry *old_head;
    do {
        old_head = hash_table[key];
        newent->next = old_head;
    } while (!atomic_compare_exchange_weak(&hash_table[key], &old_head, newent));
}

bool ht_search(int i, int j, int matrix, int *score) {
    unsigned int key = hash(i, j, matrix);
    ht_entry *entry = atomic_load(&hash_table[key]);
    while (entry) {
        if (keymatch(entry, i, j, matrix)) {
            *score = entry->score;
            return true;
        }
        entry = entry->next;
    }
    return false;
}

// Cleanup function
void cleanup_hash_table() {
    for (int i = 0; i < TABLE_SIZE; i++) {
        ht_entry *entry = hash_table[i];
        while (entry) {
            ht_entry *next = entry->next;
            free(entry);
            entry = next;
        }
    }
}