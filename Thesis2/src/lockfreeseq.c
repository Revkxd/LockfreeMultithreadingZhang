#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

#include "constants.h"
#include "tables.h"
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

typedef struct t_data {
    int thread_id;
    /*
        (i, j) for thread to start, no need for matrix because the algorithm flow computes for it
    */
    int i;
    int j;
    int score;
    char* dnaSequence;
    char* proteinSequence;
    int N;
    int M;
    int gep;
    int gop;
    int frameshift_penalty;
} thread_data_t;

#define MAX_THREADS 4
#define MASTER_THREAD_ID 0
static unsigned int defined_max_threads = MAX_THREADS;
static thread_data_t thread_data[MAX_THREADS];
static pthread_t threads[MAX_THREADS];

static pthread_mutex_t term_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t term_cond = PTHREAD_COND_INITIALIZER;
static int term_thread_id = -1;

static unsigned int num_active_threads = 0;
int active_thread_ids[MAX_THREADS];

int calculateI(char* dnaSequence, char* proteinSequence, int i, int j, int gep, int gop, int frameshift_penalty);
int calculateD(char* dnaSequence, char* proteinSequence, int i, int j, int gep, int gop, int frameshift_penalty);
int calculateC(char* dnaSequence, char* proteinSequence, int i, int j, int gep, int gop, int frameshift_penalty);

// void initializations(char* dna, char* protein, int N, int M, int gep, int gop, int frameshift_penalty) {
//     return;
//     int j;
//     for(j = 0; j < M + 1; j++) {
//         ht_insert(0, j, 1, -999);
//         ht_insert(0, j, 2, -999);
//         ht_insert(2, j, 2, -999);
//         ht_insert(3, j, 2, -999);
//         ht_insert(1, j, 2, calculateC(dna, protein, 0, j, gep, gop, frameshift_penalty) - gop - gep);
//     }

//     for(j = 1; j < M + 1; j++) {
//         ht_insert(0, j, 3, 0);
//         ht_insert(j, 0, 3, 0);
//         ht_insert(1, j, 3, max_of_three(calculateI(dna, protein, 1, j, gep, gop, frameshift_penalty),
//                                         calculateD(dna, protein, 1, j, gep, gop, frameshift_penalty),
//                                         calculateC(dna, protein, 0, j - 1, gep, gop, frameshift_penalty) + get_score(protein[j - 1], get_translated_codon(dna, 1))));
//         ht_insert(2, j, 3, max_of_two(calculateI(dna, protein, 2, j, gep, gop, frameshift_penalty),
//                                     calculateC(dna, protein, 0, j - 1, gep, gop, frameshift_penalty) + get_score(protein[j - 1], get_translated_codon(dna, 2)) - frameshift_penalty));
//         ht_insert(3, j, 3, max_of_two(calculateI(dna, protein, 3, j, gep, gop, frameshift_penalty),
//                                     calculateC(dna, protein, 1, j - 1, gep, gop, frameshift_penalty) + get_score(protein[j - 1], get_translated_codon(dna, 3)) - frameshift_penalty));
//     }
// }

int calculateI(char* dnaSequence, char* proteinSequence, int i, int j, int gep, int gop, int frameshift_penalty) {
    int score = -999;
    if (ht_search(i, j, 1, &score)) {
        return score;
    }
    // Base case: if j is 0, return 0
    if (j == 0) {
        return -999;
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
        if (i != 1) {
            return -999;
        }
        else {
            return calculateC(dnaSequence, proteinSequence, 0, j, gep, gop, frameshift_penalty) - gop - gep;
        }
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
        ht_insert(i, j, 3, 0);
        return 0;
    }
    else if (i == 1 && j > 0) {
        score = max_of_three(
            calculateI(dnaSequence, proteinSequence, 1, j, gep, gop, frameshift_penalty),
            calculateD(dnaSequence, proteinSequence, 1, j, gep, gop, frameshift_penalty),
            calculateC(dnaSequence, proteinSequence, 0, j-1, gep, gop, frameshift_penalty) + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 1))
        );
        score = score < 0 ? 0 : score;
        ht_insert(i, j, 3, score);
        return score;
    }
    else if (i == 2 && j > 0) {
        score = max_of_two(
            calculateI(dnaSequence, proteinSequence, 2, j, gep, gop, frameshift_penalty),
            calculateC(dnaSequence, proteinSequence, 0, j-1, gep, gop, frameshift_penalty) + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 2)) - frameshift_penalty
        );
        score = score < 0 ? 0 : score;
        ht_insert(i, j, 3, score);
        return score;
    }
    else if (i == 3 && j > 0) {
        score = max_of_two(
            calculateI(dnaSequence, proteinSequence, 3, j, gep, gop, frameshift_penalty),
            calculateC(dnaSequence, proteinSequence, 1, j-1, gep, gop, frameshift_penalty) + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 3)) - frameshift_penalty
        );
        score = score < 0 ? 0 : score;
        ht_insert(i, j, 3, score);
        return score;
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

int modded_three_frame(char* dnaSequence, char* proteinSequence, int start_I, int start_J, int N, int M, int gep, int gop, int frameshift_penalty) {
    int i, j;
    int max_val = -999;
    int C[N][M+1];

    max_val = calculateC(dnaSequence, proteinSequence, start_I, start_J, gep, gop, frameshift_penalty);

    // Matrix Filling
    // for(i = N-1; i >= 4; i--) {
        // for(j = M; j > 0; j--) {
    for(i = 4; i < N-1; i++) {
        for(j = 1; j < M+1; j++) {
            C[i][j] = calculateC(dnaSequence, proteinSequence, i, j, gep, gop, frameshift_penalty);
        }
    }

    for(i = 0; i < N; i++) {
        for(j = 0; j < M + 1; j++) {
            max_val = C[i][j] > max_val ? C[i][j] : max_val;
        }
    }
    
    return max_val;
}

void *three_frame_thread(void *arg) {
    thread_data_t *data = (thread_data_t *)arg;
    int i = data->i;
    int j = data->j;
    int N = data->N;
    int M = data->M;
    int I[N][M + 1], D[N][M + 1], C[N][M+1];
    int TI[N][M + 1], TD[N][M + 1], TC[N][M + 1];
    
    pthread_mutex_lock(&term_mutex);
    data->score = modded_three_frame(data->dnaSequence, data->proteinSequence, data->i, data->j, data->N, data->M, data->gep, data->gop, data->frameshift_penalty);
    if(term_thread_id == -1) {
        term_thread_id = data->thread_id;
    }
    pthread_cond_signal(&term_cond);
    pthread_mutex_unlock(&term_mutex);
    return &data->score;
}

static int three_frame_master_thread(char* dnaSequence, char* proteinSequence, int N, int M, int gep, int gop, int frameshift_penalty) {
    int actindex=0;
    int new_thread_id;
    int finished_thread_id;
    int rc;
    int *score;

    term_thread_id = -1;
    num_active_threads = 0;

    srand(time(NULL));
    while (num_active_threads < defined_max_threads)
    {
        new_thread_id = num_active_threads;
        thread_data[num_active_threads].thread_id = num_active_threads;
        thread_data[num_active_threads].i = rand() % N;
        thread_data[num_active_threads].j = rand() % M;
        #ifdef PRINTERS
        printf("i, j for thread %d: %d, %d\n", new_thread_id, thread_data[num_active_threads].i, thread_data[num_active_threads].j);
        #endif
        thread_data[num_active_threads].score = -999;
        thread_data[num_active_threads].dnaSequence=dnaSequence;
        thread_data[num_active_threads].proteinSequence=proteinSequence;
        thread_data[num_active_threads].N = N;
        thread_data[num_active_threads].M = M;
        thread_data[num_active_threads].gep = gep;
        thread_data[num_active_threads].gop = gop;
        thread_data[num_active_threads].frameshift_penalty = frameshift_penalty;
        rc = pthread_create(&threads[num_active_threads], NULL, three_frame_thread, (void *)&thread_data[num_active_threads]);
        if (rc)
        {
            printf("ERROR; return code from pthread_create() for thread_id %d is %d\n", new_thread_id, rc);
            exit(-1);
        }
        active_thread_ids[actindex++] = num_active_threads;
        num_active_threads++;
    }
    pthread_mutex_lock(&term_mutex);
    if (term_thread_id == -1)
    {
        pthread_cond_wait(&term_cond, &term_mutex);
    }
    finished_thread_id = term_thread_id;
    pthread_mutex_unlock(&term_mutex);
    rc = pthread_join(threads[finished_thread_id], (void **)&score);
    if (rc)
    {
        printf("ERROR; return code from pthread_join() is %d\n", rc);
        exit(-1);
    }
    return *score;
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
    init_hash_table();
    // initializations(dnaSequence, proteinSequence, N, M, gep, gop, frameshift_penalty);
    // max1 = modded_three_frame(dnaSequence, proteinSequence, N, M, I, D, C, gep, gop, frameshift_penalty);
    max1 = three_frame_master_thread(dnaSequence, proteinSequence, N, M, gep, gop, frameshift_penalty);

    #ifdef PRINTERS
    printf("Reverse Complement:\n");
    #endif
    init_hash_table();
    // initializations(dnaSequence, proteinSequence, N, M, gep, gop, frameshift_penalty);
    reverse_complement(dnaSequence);
    max2 = three_frame_master_thread(dnaSequence, proteinSequence, N, M, gep, gop, frameshift_penalty);

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
    #ifndef ITER
    init_hash_table();
    #endif
    printf("DNA Sequence: %s\n", dna);
    printf("Protein Sequence: %s\n", protein);
    start = clock();
    printf("Score: %d\n\n", six_frame(dna, protein));
    end = clock();
    time_taken = (double)(end - start)*1e3 / CLOCKS_PER_SEC;
    printf("Run time taken: %f ms\n\n", time_taken);

    return 0;
}