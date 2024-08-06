#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdatomic.h>

#define TABLE_SIZE 500000

typedef struct ht_ent {
    struct ht_ent *next;
    int i;
    int j;
    int score;
} ht_entry;

static ht_entry *I_table[TABLE_SIZE];
static ht_entry *D_table[TABLE_SIZE];
static ht_entry *C_table[TABLE_SIZE];

void init_hash_table() {
    for (int i = 0; i < TABLE_SIZE; i++) {
        I_table[i] = NULL;
        D_table[i] = NULL;
        C_table[i] = NULL;
    }
}

static unsigned long hash6432shift(unsigned long long key)
{
  key = (~key) + (key << 18); /* key = (key << 18) - key - 1; */
  key = key ^ (key >> 31);
  key = key * 21; /* key = (key + (key << 2)) + (key << 4); */
  key = key ^ (key >> 11);
  key = key + (key << 6);
  key = key ^ (key >> 22);
  return (unsigned long) key;
}

unsigned long hash(int i, int j) {
    unsigned long long key; //= (i << 16) | (j & 0xffff);
    key = ((i << 16) | j) & 0xFFFFFFFF;
    key = hash6432shift(key);
    return key % TABLE_SIZE;
}

int keymatch(ht_entry *entry, int i, int j) {
    return entry->i == i && entry->j == j;
}

void ht_insert(int i, int j, int matrix, int score) {
    unsigned long key = hash(i, j);
    ht_entry *newent = (ht_entry *)malloc(sizeof(ht_entry));
    if (!newent) {
        fprintf(stderr, "Memory allocation failed\n");
        return;
    }

    newent->i = i;
    newent->j = j;
    newent->score = score;
    newent->next = NULL;

    ht_entry *old_head;
    switch(matrix) {
        case 1:
            do {
                old_head = I_table[key];
                newent->next = old_head;
            } while (!atomic_compare_exchange_weak(&I_table[key], &old_head, newent));
            break;
        case 2:
            do {
                old_head = D_table[key];
                newent->next = old_head;
            } while (!atomic_compare_exchange_weak(&D_table[key], &old_head, newent));
            break;
        case 3:
            do {
                old_head = C_table[key];
                newent->next = old_head;
            } while (!atomic_compare_exchange_weak(&C_table[key], &old_head, newent));
            break;
    }
}

bool ht_search(int i, int j, int matrix, int *score) {
    unsigned long key = hash(i, j);
    ht_entry *entry;
    switch(matrix) {
        case 1: entry = atomic_load(&I_table[key]); break;
        case 2: entry = atomic_load(&D_table[key]); break;
        case 3: entry = atomic_load(&C_table[key]);
    }
    while (entry) {
        if (keymatch(entry, i, j)) {
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
        ht_entry *I_entry = I_table[i];
        while (I_entry) {
            ht_entry *next = I_entry->next;
            free(I_entry);
            I_entry = next;
        }

        ht_entry *D_entry = D_table[i];
        while (D_entry) {
            ht_entry *next = D_entry->next;
            free(D_entry);
            D_entry = next;
        }

        ht_entry *C_entry = C_table[i];
        while (C_entry) {
            ht_entry *next = C_entry->next;
            free(C_entry);
            C_entry = next;
        }
    }
}