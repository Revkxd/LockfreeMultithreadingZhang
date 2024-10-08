#define OWNHT
#ifdef OWNHT
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

#include "../src/constants.h"
#include "../src/converters.c"
#include "../src/maxfunctions.c"
#include "selfht2.c"

#define NUM_INSTANCES 10
#define PRINTERS 1
#undef PRINTERS
#define THREADING 1

int reversed = 0;
entry_t *orig_ht[TABLE_SIZE];
entry_t *reverse_ht[TABLE_SIZE];

typedef struct {
    char *dnaSequence;
    char *proteinSequence;
    entry_t **hash_table;
    int *max_val;
} thread_data;

typedef struct {
    int thread_no;
    char *dnaSequence;
    char *proteinSequence;
    int best_score;
    double start;
    double end;
} main_thread_data;

void modded_three_frame(char* dnaSequence, char* proteinSequence, entry_t *hash_table[], int* max_value);

void* middleman(void *args) {
    thread_data *data = (thread_data *)args;
    modded_three_frame(data->dnaSequence, data->proteinSequence, data->hash_table, data->max_val);
}

void modded_three_frame(char* dnaSequence, char* proteinSequence, entry_t *hash_table[], int* max_value) {
    int N = strlen(dnaSequence), M = strlen(proteinSequence);
    int gep = 2, gop = 3, frameshift_penalty = 4;
    int I[N][M + 1], D[N][M + 1], C[N][M+1];
    int TI[N][M + 1], TD[N][M + 1], TC[N][M + 1];
    int i,j;

    // Anti-Garbage Values
    for(i = 0; i < N; i++) {
        for(j = 0; j < M + 1; j++) {
            I[i][j] = 0;
            D[i][j] = 0;
            C[i][j] = 0;
            TI[i][j] = 0;
            TD[i][j] = 0;
            TC[i][j] = 0;
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

    #ifdef FILE_TRACE
    char filename[50];
    char letter;
    #endif
    // Matrix Filling
    for(i = 0; i < N; i++) {
        #ifdef FILE_TRACE
        letter = 'A' + i;
        strcpy(filename, &letter);
        if (reversed)
            strcat(filename, &letter);
        strcat(filename, ".txt");
        FILE *fp = fopen(filename, "w");
        #endif
        for(j = 1; j < M + 1; j++) {
            if(!ht_search(i, j, 1, &I[i][j], hash_table)) {
                I[i][j] = max_of_two(I[i][j-1] - gep, C[i][j-1] - gop - gep);
                ht_insert(i, j, 1, I[i][j], hash_table);
            }
            if (i < 4)
                continue;
            if(!ht_search(i, j, 2, &D[i][j], hash_table)) {
                D[i][j] = max_of_two(D[i-3][j] - gep, C[i-3][j] - gop - gep);
                ht_insert(i, j, 2, D[i][j], hash_table);
            }
            if(!ht_search(i, j, 3, &C[i][j], hash_table)) {
                C[i][j] = max_of_five(
                    I[i][j],
                    D[i][j],
                    C[i-2][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, i)) - frameshift_penalty,
                    C[i-3][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, i)),
                    C[i-4][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, i)) - frameshift_penalty
                );
                ht_insert(i, j, 3, C[i][j], hash_table);
            }
        }
        #ifdef FILE_TRACE
        print_table_to_file(fp, hash_table);
        fclose(fp);
        #endif
    }

    // Traceback Matrix Filling
    for(i = 0; i < N; i++) {
        #ifdef FILE_TRACE
        letter = 'A' + i;
        strcpy(filename, &letter);
        if (reversed)
            strcat(filename, &letter);
        strcat(filename, "_traceback.txt");
        FILE *fp = fopen(filename, "w");
        #endif
        for(j = 0; j < M; j++) {
            if(!ht_search(i, j, 4, &TI[i][j], hash_table)) {
                TI[i][j] = I[i][j] == I[i][j-1] - gep ? 0 : (I[i][j] == C[i][j-1] - gop - gep ? 1 : -999);
                ht_insert(i, j, 4, TI[i][j], hash_table);
            }
            if(!ht_search(i, j, 5, &TD[i][j], hash_table)) {
                TD[i][j] = D[i][j] == D[i-3][j] - gep ? 0 : (D[i][j] == C[i-3][j] - gop - gep ? 1 : -999);
                ht_insert(i, j, 5, TD[i][j], hash_table);
            }
            if(!ht_search(i, j, 6, &TC[i][j], hash_table)) {
                if(C[i][j] == I[i][j])
                    TC[i][j] = -2;
                else if(C[i][j] == D[i][j])
                    TC[i][j] = -1;
                else if(C[i][j] == C[i-2][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, i)) - frameshift_penalty)
                    TC[i][j] = 2;
                else if(C[i][j] == C[i-3][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, i)))
                    TC[i][j] = 3;
                else if(C[i][j] == C[i-4][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, i)) - frameshift_penalty)
                    TC[i][j] = 4;
                else
                    TC[i][j] = 0;
                ht_insert(i, j, 6, TC[i][j], hash_table);
            }
        }
        #ifdef FILE_TRACE
        print_table_to_file(fp, hash_table);
        fclose(fp);
        #endif
    }

    // Print the matrices for debugging
    #ifdef PRINTERS
    printf("I Matrix:\n");
    print_matrix(N, M + 1, I);
    printf("D Matrix:\n");
    print_matrix(N, M + 1, D);
    printf("C Matrix:\n");
    print_matrix(N, M + 1, C);
    printf("TI Matrix:\n");
    print_matrix(N, M + 1, TI);
    printf("TD Matrix:\n");
    print_matrix(N, M + 1, TD);
    printf("TC Matrix:\n");
    print_matrix(N, M + 1, TC);
    #endif

    int max_val = -999;
    for(i = 0; i < N; i++) {
        for(j = 0; j < M + 1; j++) {
            max_val = C[i][j] > max_val ? C[i][j] : max_val;
        }
    }
    *max_value = max_val;   
    // return max_val;
}

int six_frame(char* dnaSequence, char* proteinSequence);
void* six_frame_middleman(void* args) {
    main_thread_data *data = (main_thread_data *)args;
    data->start = clock();
    data->best_score = six_frame(data->dnaSequence, data->proteinSequence);
    data->end = clock();
}

int six_frame(char* dnaSequence, char* proteinSequence) {
    int N = strlen(dnaSequence), M = strlen(proteinSequence);
    int C1[N][M + 1], C2[N][M + 1];
    int max1, max2;
    pthread_t thread1, thread2;
    // int retval[2];

    #ifdef PRINTERS
    printf("First Run:\n");
    #endif
    reversed = 0;
    thread_data *data1 = (thread_data *)malloc(sizeof(thread_data));
    data1->dnaSequence = dnaSequence;
    data1->proteinSequence = proteinSequence;
    data1->hash_table = orig_ht;
    data1->max_val = &max1;
    pthread_create(&thread1, NULL, middleman, (void *)data1);
    // modded_three_frame(dnaSequence, proteinSequence, orig_ht, &retval[0]);

    #ifdef PRINTERS
    printf("Reverse Complement:\n");
    #endif
    reversed = 1;
    char *dnaSequence2 = (char *)malloc(N * sizeof(char));
    strcpy(dnaSequence2, dnaSequence);
    reverse_complement(dnaSequence2);
    thread_data *data2 = (thread_data *)malloc(sizeof(thread_data));
    data2->dnaSequence = dnaSequence2;
    data2->proteinSequence = proteinSequence;
    data2->hash_table = reverse_ht;
    data2->max_val = &max2;
    pthread_create(&thread2, NULL, middleman, (void *)data2);
    // modded_three_frame(dnaSequence2, proteinSequence, reverse_ht, &retval[1]);

    pthread_join(thread1, NULL);
    pthread_join(thread2, NULL);
    // max1 = retval[0];
    // max2 = retval[1];

    // printf("Value of max1, max2: %d, %d\n", max1, max2);
    return max_of_two(max1, max2);
}

int main() {
    String dnaSequences[] = {"ATTGACAACCGCGTCCGCCGCCGCTTCAAGGGCCAGTACTTGATGCCCAACATTGGCTACGGCTCCAACAAGCGCACCCGCCACATGTTGCCCACCGGCT", "GCTACGTCCGCTCCTCCATGTCCTTGTCCGGCTACATGCCCCCCTTGTGCGACCCCAAGGACGGCCACTTGTTGTTGGACGGCGGCTACGTCAACAACT", "GAGCCCACCTCCGAGATTTTGCAGAACCCCGCCCGCGTCTTGCGCCAGCAGTTGAAGGTCTTGTCCGTCATTGACGGCCAGTCCTACGAGCCCTTGAAGG", "CCCGGCGCCGGCTCCGGCCACGGCCACGGCCCCAACGGCGGCTCCAACTCCTCCTCCTGCACCCCCCCCTCCTCCAACCCCCACATTACCGGCTACGTCG"};
    String proteinSequences[] = {"IDNRVRRRFKGQYLMPNIGYGSNKRTRHMLPTGF", "RYVRSSMSLSGYMPPLCDPKDGHLLLDGGYVNNL", "EPTSEILQNPARVLRQQLKVLSVIDGQSYEPLKD", "PGAGSGHGHGPNGGSNSSSCTPPSSNPHITGYVD"};
    int i;
    double time_taken, start, end;
    #ifdef THREADING
    main_thread_data *data[NUM_INSTANCES];
    for(i = 0; i < NUM_INSTANCES; i++) {
        data[i] = (main_thread_data *)malloc(sizeof(main_thread_data));
    }
    pthread_t threads[NUM_INSTANCES];
    #endif
    for(i = 0; i < 4; i++) {
        init_hash_table(orig_ht);
        init_hash_table(reverse_ht);
        printf("DNA Sequence: %s\n", dnaSequences[i]);
        printf("Protein Sequence: %s\n", proteinSequences[i]);
        start = clock();
        #ifndef THREADING
        printf("Score: %d\n\n", six_frame(dnaSequences[i], proteinSequences[i]));
        #else
        for(int j = 0; j < NUM_INSTANCES; j++) {
            data[j]->thread_no = j;
            data[j]->dnaSequence = dnaSequences[i];
            data[j]->proteinSequence = proteinSequences[i];
            pthread_create(&threads[j], NULL, six_frame_middleman, (void *)data[j]);
        }
        for(int j = 0; j < NUM_INSTANCES; j++) {
            pthread_join(threads[j], NULL);
            printf("Thread %d: %d\n", j, data[j]->best_score);
            printf("\tTime taken: %f ms\n", (data[j]->end - data[j]->start)*1e3 / CLOCKS_PER_SEC);
        }
        #endif
        end = clock();
        time_taken = (double)(end - start)*1e3 / CLOCKS_PER_SEC;
        printf("Run %d time taken: %f ms\n\n", i, time_taken);
        // break;
    }
    return 0;
}
#endif