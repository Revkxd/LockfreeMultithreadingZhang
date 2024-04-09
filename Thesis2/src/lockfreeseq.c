#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

#include "maxfunctions.c"
#include "converters.c"
#include "splitlfht.c"

#define PRINTERS 1
#undef PRINTERS

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

void initializations(char* dna, char* protein, int N, int M, int gep, int gop, int frameshift_penalty) {
    return; // disabled because for some reason when done the answers are wrong.
    int j;
    for(j = 0; j < M + 1; j++) {
        ht_insert(0, j, 1, -999);
        ht_insert(0, j, 2, -999);
        ht_insert(2, j, 2, -999);
        ht_insert(3, j, 2, -999);
        ht_insert(1, j, 2, calculateC(dna, protein, 0, j, gep, gop, frameshift_penalty) - gop - gep);
    }

    for(j = 1; j < M + 1; j++) {
        ht_insert(0, j, 3, 0);
        ht_insert(j, 0, 3, 0);
        ht_insert(1, j, 3, max_of_three(calculateI(dna, protein, 1, j, gep, gop, frameshift_penalty),
                                        calculateD(dna, protein, 1, j, gep, gop, frameshift_penalty),
                                        calculateC(dna, protein, 0, j - 1, gep, gop, frameshift_penalty) + get_score(protein[j - 1], get_translated_codon(dna, 1))));
        ht_insert(2, j, 3, max_of_two(calculateI(dna, protein, 2, j, gep, gop, frameshift_penalty),
                                    calculateC(dna, protein, 0, j - 1, gep, gop, frameshift_penalty) + get_score(protein[j - 1], get_translated_codon(dna, 2)) - frameshift_penalty));
        ht_insert(3, j, 3, max_of_two(calculateI(dna, protein, 3, j, gep, gop, frameshift_penalty),
                                    calculateC(dna, protein, 1, j - 1, gep, gop, frameshift_penalty) + get_score(protein[j - 1], get_translated_codon(dna, 3)) - frameshift_penalty));
        ht_insert(4, j, 3, max_of_four(calculateI(dna, protein, 4, j, gep, gop, frameshift_penalty),
                                    calculateD(dna, protein, 4, j, gep, gop, frameshift_penalty),
                                    calculateC(dna, protein, 1, j - 1, gep, gop, frameshift_penalty) + get_score(protein[j - 1], get_translated_codon(dna, 4)),
                                    calculateC(dna, protein, 2, j - 1, gep, gop, frameshift_penalty) + get_score(protein[j - 1], get_translated_codon(dna, 4)) - frameshift_penalty));
    }
}

void matrixInitialize(char* dnaSequence, char* proteinSequence, int N, int M, int I[][M+1], int D[][M+1], int C[][M+1], int TI[][M+1], int TD[][M+1], int TC[][M+1], int gep, int gop, int frameshift_penalty) {
    return;
    int i, j;
    // Initialization
    for(j = 0; j < M + 1; j++) {
        I[0][j] = -999;
        D[0][j] = -999;
        D[2][j] = -999;
        D[3][j] = -999;
        D[1][j] = C[0][j] - gop - gep;
        TI[j][0] = -999;
        TD[0][j] = -999;
        TD[2][j] = -999;
        TD[3][j] = -999;
        TD[1][j] = 1;
    }

    C[0][0] = 0;
    TC[0][0] = 0;
    for(j = 1; j < M + 1; j++) {
        C[0][j] = 0;
        C[j][0] = 0;
        C[1][j] = max_of_three(I[1][j], D[1][j], C[0][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 1)));
        C[2][j] = max_of_two(I[2][j], C[0][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 2)) - frameshift_penalty);
        C[3][j] = max_of_two(I[3][j], C[1][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 3)) - frameshift_penalty);
        C[4][j] = max_of_four(I[4][j], D[4][j], C[1][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 4)), C[2][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 4)) - frameshift_penalty);

        TC[0][j] = 0;
        TC[j][0] = 0;
        TC[1][j] = C[1][j] == I[1][j] ? -2 : (C[1][j] == D[1][j] ? -1 : (C[1][j] == C[0][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 1)) ? 1 : 0));
        TC[2][j] = C[2][j] == I[2][j] ? -2 : (C[2][j] == C[0][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 2)) - frameshift_penalty ? 2 : 0);
        TC[3][j] = C[3][j] == I[3][j] ? -2 : (C[3][j] == C[1][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 3)) - frameshift_penalty ? 2 : 0);
        if(C[4][j] == I[4][j])
            TC[4][j] = -2;
        else if(C[4][j] == D[4][j])
            TC[4][j] = -1;
        else if(C[4][j] == C[1][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 4)))
            TC[4][j] = 3;
        else if(C[4][j] == C[2][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 4)) - frameshift_penalty)
            TC[4][j] = 2;
        else
            TC[4][j] = 0;
    }

    // Change all negative values in C matrix to 0
    for(i = 0; i < N; i++) {
        for(j = 0; j < M + 1; j++) {
            if(C[i][j] < 0)
                C[i][j] = 0;
        }
    }
}

int calculateI(char* dnaSequence, char* proteinSequence, int i, int j, int gep, int gop, int frameshift_penalty) {
    int score = -999;
    if (ht_search(i, j, 1, &score)) {
        return score;
    }
    // Base case: if j is 0, return 0
    if (i == 0 || j == 0) {
        return 0;
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
        return 0;
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
        return 0;
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

int modded_three_frame(char* dnaSequence, char* proteinSequence, int N, int M, int I[][M+1], int D[][M+1], int C[][M+1], int TI[][M+1], int TD[][M+1], int TC[][M+1], int gep, int gop, int frameshift_penalty) {
    int i, j;

    // Matrix Filling
    for(i = 0; i < N; i++) {
        for(j = 1; j < M + 1; j++) {
            // I[i][j] = calculateI(dnaSequence, proteinSequence, i, j, gep, gop, frameshift_penalty);
            if (i < 4)
                continue;
            // D[i][j] = calculateD(dnaSequence, proteinSequence, i, j, gep, gop, frameshift_penalty);
            C[i][j] = calculateC(dnaSequence, proteinSequence, i, j, gep, gop, frameshift_penalty);
        }
    }

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

    int max_val = -999;
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
    data->score = modded_three_frame(data->dnaSequence, data->proteinSequence, data->N, data->M, I, D, C, TI, TD, TC, data->gep, data->gop, data->frameshift_penalty);
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
    initializations(dnaSequence, proteinSequence, N, M, gep, gop, frameshift_penalty);
    // max1 = modded_three_frame(dnaSequence, proteinSequence, N, M, I, D, C, gep, gop, frameshift_penalty);
    max1 = three_frame_master_thread(dnaSequence, proteinSequence, N, M, gep, gop, frameshift_penalty);

    #ifdef PRINTERS
    printf("Reverse Complement:\n");
    #endif
    init_hash_table();
    initializations(dnaSequence, proteinSequence, N, M, gep, gop, frameshift_penalty);
    reverse_complement(dnaSequence);
    // max2 = modded_three_frame(dnaSequence, proteinSequence, N, M, I, D, C, gep, gop, frameshift_penalty);
    max2 = three_frame_master_thread(dnaSequence, proteinSequence, N, M, gep, gop, frameshift_penalty);

    return max_of_two(max1, max2);
}