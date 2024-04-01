#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

#include "maxfunctions.c"
#include "converters.c"
#include "lfhtfixed.c"

#define PRINTERS 1
#undef PRINTERS

time_t overall_start;
time_t current_clock;
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

#define MAX_THREADS 6
#define MASTER_THREAD_ID 0
static unsigned int defined_max_threads = MAX_THREADS;
static thread_data_t thread_data[MAX_THREADS];
static pthread_t threads[MAX_THREADS];

static pthread_mutex_t term_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t term_cond = PTHREAD_COND_INITIALIZER;
static int term_thread_id = -1;

static unsigned int num_active_threads = 0;
int active_thread_ids[MAX_THREADS];

FILE *starts;
FILE *insert_I;
FILE *insert_D;
FILE *insert_C;

int matrixInitialize(char* dnaSequence, char* proteinSequence, int N, int M, int I[][M+1], int D[][M+1], int C[][M+1], int gep, int gop, int frameshift_penalty) {
    int i, j;
    // Anti-Garbage Values
    for(i = 0; i < N; i++) {
        for(j = 0; j < M + 1; j++) {
            I[i][j] = 0;
            D[i][j] = 0;
            C[i][j] = 0;
        }
    }

    // NOTE: INFINITY is a double, just converted to int to save memory
    // Initialization
    for(j = 0; j < M + 1; j++) {
        I[0][j] = -999;
        D[0][j] = -999;
        D[2][j] = -999;
        D[3][j] = -999;
        D[1][j] = C[0][j] - gop - gep;
    }

    C[0][0] = 0;
    for(j = 1; j < M + 1; j++) {
        C[0][j] = 0;
        C[j][0] = 0;
        C[1][j] = max_of_three(I[1][j], D[1][j], C[0][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 1)));
        C[2][j] = max_of_two(I[2][j], C[0][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 2)) - frameshift_penalty);
        C[3][j] = max_of_two(I[3][j], C[1][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 3)) - frameshift_penalty);
        C[4][j] = max_of_four(I[4][j], D[4][j], C[1][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 4)), C[2][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 4)) - frameshift_penalty);
    }

    // Change all negative values in C matrix to 0
    for(i = 0; i < N; i++) {
        for(j = 0; j < M + 1; j++) {
            if(C[i][j] < 0)
                C[i][j] = 0;
        }
    }
}

int calculateI(int thread_id, char* dnaSequence, char* proteinSequence, int i, int j, int gep, int gop, int frameshift_penalty);
int calculateD(int thread_id, char* dnaSequence, char* proteinSequence, int i, int j, int gep, int gop, int frameshift_penalty);
int calculateC(int thread_id, char* dnaSequence, char* proteinSequence, int i, int j, int gep, int gop, int frameshift_penalty);

int calculateI(int thread_id, char* dnaSequence, char* proteinSequence, int i, int j, int gep, int gop, int frameshift_penalty) {
    int score = -999;
    if (ht_search(i, j, 1, &score)) {
        fprintf(insert_I, "FOUND! Tick: %ld Thread id: %d, i: %d, j: %d, matrix: %c, score: %d\n", clock() - overall_start, thread_id, i, j, 'I', score);
        return score;
    }
    // Base case: if j is 0, return 0
    if (i == 0 || j == 0) {
        return 0;
    }

    // Recursive case for I matrix
    score = max_of_two(calculateI(thread_id, dnaSequence, proteinSequence, i, j - 1, gep, gop, frameshift_penalty) - gep, calculateC(thread_id, dnaSequence, proteinSequence, i, j - 1, gep, gop, frameshift_penalty) - gop - gep);
    fprintf(insert_I, "INSERT Tick: %ld || Thread id: %d, i: %d, j: %d, matrix: I, score: %d\n", clock() - overall_start, thread_id, i, j, score);
    ht_insert(i, j, 1, score);
    return score;
}

int calculateD(int thread_id, char* dnaSequence, char* proteinSequence, int i, int j, int gep, int gop, int frameshift_penalty) {
    int score = -999;
    if(ht_search(i, j, 2, &score)) {
        fprintf(insert_D, "FOUND! Tick: %ld Thread id: %d, i: %d, j: %d, matrix: %c, score: %d\n", clock() - overall_start, thread_id, i, j, 'D', score);
        return score;
    }
    // Base case: if i is less than 4 or j is 0, return 0
    if (i < 4 || j == 0) {
        return 0;
    }

    // Recursive case for D matrix
    score = max_of_two(calculateD(thread_id, dnaSequence, proteinSequence, i - 3, j, gep, gop, frameshift_penalty) - gep, calculateC(thread_id, dnaSequence, proteinSequence, i - 3, j, gep, gop, frameshift_penalty) - gop - gep);
    fprintf(insert_D, "INSERT Tick: %ld || Thread id: %d, i: %d, j: %d, matrix: D, score: %d\n", clock() - overall_start, thread_id, i, j, score);
    ht_insert(i, j, 2, score);
    return score;
}

int calculateC(int thread_id, char* dnaSequence, char* proteinSequence, int i, int j, int gep, int gop, int frameshift_penalty) {
    int score = -999;
    if (ht_search(i, j, 3, &score)) {
        fprintf(insert_C, "FOUND! Tick: %ld Thread id: %d, i: %d, j: %d, matrix: %c, score: %d\n", clock() - overall_start, thread_id, i, j, 'C', score);
        return score;
    }
    // Base case: if i or j is 0, return 0
    if (i == 0 || j == 0) {
        return 0;
    }

    // Recursive cases for C matrix
    int I_val = calculateI(thread_id, dnaSequence, proteinSequence, i, j, gep, gop, frameshift_penalty);
    int D_val = calculateD(thread_id, dnaSequence, proteinSequence, i, j, gep, gop, frameshift_penalty);

    // Calculate C[i][j] based on I[i][j], D[i][j], and C[i][j]
    // Adjust this part based on your actual calculation logic
    score = max_of_five(
        I_val,
        D_val,
        calculateC(thread_id, dnaSequence, proteinSequence, i - 2, j - 1, gep, gop, frameshift_penalty) + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, i)) - frameshift_penalty,
        calculateC(thread_id, dnaSequence, proteinSequence, i - 3, j - 1, gep, gop, frameshift_penalty) + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, i)),
        calculateC(thread_id, dnaSequence, proteinSequence, i - 4, j - 1, gep, gop, frameshift_penalty) + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, i)) - frameshift_penalty
    );
    fprintf(insert_C, "INSERT Tick: %ld || Thread id: %d, i: %d, j: %d, matrix: C, score: %d\n", clock() - overall_start, thread_id, i, j, score);
    ht_insert(i, j, 3, score);
    return score;
}

int modded_three_frame(int thread_id, char* dnaSequence, char* proteinSequence, int N, int M, int I[][M+1], int D[][M+1], int C[][M+1], int gep, int gop, int frameshift_penalty) {
    int i, j;

    matrixInitialize(dnaSequence, proteinSequence, N, M, I, D, C, gep, gop, frameshift_penalty);

    // Matrix Filling
    for(i = 0; i < N; i++) {
        for(j = 1; j < M + 1; j++) {
            I[i][j] = calculateI(thread_id, dnaSequence, proteinSequence, i, j, gep, gop, frameshift_penalty);
            if (i < 4)
                continue;
            D[i][j] = calculateD(thread_id, dnaSequence, proteinSequence, i, j, gep, gop, frameshift_penalty);
            C[i][j] = calculateC(thread_id, dnaSequence, proteinSequence, i, j, gep, gop, frameshift_penalty);
        }
    }

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
    matrixInitialize(data->dnaSequence, data->proteinSequence, N, M, I, D, C, data->gep, data->gop, data->frameshift_penalty);
    
    pthread_mutex_lock(&term_mutex);
    data->score = modded_three_frame(data->thread_id, data->dnaSequence, data->proteinSequence, data->N, data->M, I, D, C, data->gep, data->gop, data->frameshift_penalty);
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
        fprintf(starts, "i, j for thread %d: %d, %d\n", new_thread_id, thread_data[num_active_threads].i, thread_data[num_active_threads].j);
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
    // max1 = modded_three_frame(dnaSequence, proteinSequence, N, M, I, D, C, gep, gop, frameshift_penalty);
    overall_start = clock();
    max1 = three_frame_master_thread(dnaSequence, proteinSequence, N, M, gep, gop, frameshift_penalty);
    cleanup_hash_table();

    #ifdef PRINTERS
    printf("Reverse Complement:\n");
    #endif
    init_hash_table();
    reverse_complement(dnaSequence);
    overall_start = clock();
    // max2 = modded_three_frame(dnaSequence, proteinSequence, N, M, I, D, C, gep, gop, frameshift_penalty);
    max2 = three_frame_master_thread(dnaSequence, proteinSequence, N, M, gep, gop, frameshift_penalty);
    cleanup_hash_table();

    return max_of_two(max1, max2);
}

int main() {
    String dnaSequences[] = {"ATTGACAACCGCGTCCGCCGC","ATTGACAACCGCGTCCGCCGCCGCTTCAAGGGCCAGTACTTGATGCCCAACATTGGCTACGGCTCCAACAAGCGCACCCGCCACATGTTGCCCACCGGCT", "GCTACGTCCGCTCCTCCATGTCCTTGTCCGGCTACATGCCCCCCTTGTGCGACCCCAAGGACGGCCACTTGTTGTTGGACGGCGGCTACGTCAACAACT", "GAGCCCACCTCCGAGATTTTGCAGAACCCCGCCCGCGTCTTGCGCCAGCAGTTGAAGGTCTTGTCCGTCATTGACGGCCAGTCCTACGAGCCCTTGAAGG", "CCCGGCGCCGGCTCCGGCCACGGCCACGGCCCCAACGGCGGCTCCAACTCCTCCTCCTGCACCCCCCCCTCCTCCAACCCCCACATTACCGGCTACGTCG"};
    String proteinSequences[] = {"IDNRVR","IDNRVRRRFKGQYLMPNIGYGSNKRTRHMLPTGF", "RYVRSSMSLSGYMPPLCDPKDGHLLLDGGYVNNL", "EPTSEILQNPARVLRQQLKVLSVIDGQSYEPLKD", "PGAGSGHGHGPNGGSNSSSCTPPSSNPHITGYVD"};
    int i;
    double time_taken, start, end;
    starts = fopen("starts.txt", "w");
    insert_I = fopen("insert_I.txt", "w");
    insert_D = fopen("insert_D.txt", "w");
    insert_C = fopen("insert_C.txt", "w");
    for(i = 0; i < 5; i++) {
        // init_hash_table();
        printf("DNA Sequence: %s\n", dnaSequences[i]);
        printf("Protein Sequence: %s\n", proteinSequences[i]);
        printf("N: %ld, M: %ld\n", strlen(dnaSequences[i]), strlen(proteinSequences[i]));
        start = clock();
        printf("Score: %d\n\n", six_frame(dnaSequences[i], proteinSequences[i]));
        end = clock();
        time_taken = (double)(end - start)*1e3 / CLOCKS_PER_SEC;
        printf("Run %d time taken: %f ms\n\n", i, time_taken);
        break;
    }
    fclose(starts);
    fclose(insert_I);
    fclose(insert_D);
    fclose(insert_C);

    // String dnaSeq = "GGCGTGGCGCAGGCGCAGAGAGGCGCACCGCGCCGGCGCAGGCGCAGAGACACATGCTAGCGCGTCCAGGGGTGGAGGCGTGGCGCAGGCGCAGAGACGCAAGCCTACGGGCGGGGGTTGGGGGGGCGTGTGTTGCAGGAGCAAAGTCGCACGGCGCCGGGCTGGGGCGGGGGGAGGGTGGCGCCGTGCACGCGCAGAAACTCACGTCACGGTGGCGCGGCGCAGAGACGGGTAGAACCTCAGTAATCCGAAAAGCCGGGATCGACCGCCCCTTGCTTGCAGCCGGGCACTACAGGACCC";
    // String protSeq = "PROHISARGVALARGVALSERPROARGGLYALAALAALASERALASERLEUCYSTHRILEALAGLNVALPROTHRSERALAPRORGLYVALARGMETPROALAPRONPROALAHISASNVALLEUVALSERALACYSARGGLYPROTHRPROPROPROSERHISARGGLYTHRCYSALASERLEUSERALAVAPRORARGARGVALSERALAHILEUGLYVALILEARGLEUPHEGLYPROSERTRPARGGLYTHRASNVALGLYPROCYSPROGLY";
    
    // char dnaSeq[STRING_MAX];
    // strcpy(dnaSeq,
    //       "ATTGACAACCGCGTCCGCCGC"
    //       );
    // char protSeq[STRING_MAX];
    // strcpy(protSeq,
    //       "IDNRVR"
    //       );
    // start = clock();
    // printf("Score: %d\n\n", six_frame(dnaSeq, protSeq));
    // end = clock();
    // time_taken = (double)(end - start)*1e3 / CLOCKS_PER_SEC;
    // printf("Run %d time taken: %f ms\n\n", i, time_taken);
    return 0;
}