#include "converters.c"
#include "../Lockfree/ht.h"

#define PRINTERS 1
// #undef PRINTERS

typedef struct tuple3_s
{
  unsigned int i;
  unsigned int j;
  char matrix;
} tuple3_t;

/*
 * ht_hash()
 *
 * Compute hash value for the table.
 * We will simply shove the low 16 bits of each index in the tuple
 * into a word and return the value modulo the tablesize.
 *
 * Paramters:
 *     key - ptr to 2-tuple to compute hash value for
 *
 * Return value:
 *     hash value
 */
static unsigned int ht_hash(const void *vkey) 
{
  int hashval;
  const tuple3_t *key = (const tuple3_t *)vkey;

  hashval = (key->i << 16) | (key->j & 0xffff);
  return hashval % HT_SIZE; 

}

/*
 * ht_keymatch()
 *
 * Compare two key structures
 *
 * Parameters:
 *    s1 - key struct ptr for first
 *    s2 - key struct ptr for second
 *
 * Return value:
 *    nonzero if s1 and s2 are equal, else 0.
 *
 */
static int ht_keymatch(const void *vs1, const void *vs2)
{
  const tuple3_t *s1 = (const tuple3_t *)vs1;
  const tuple3_t *s2 = (const tuple3_t *)vs2;

  return (s1->i == s2->i && s1->j == s2->j && s1->matrix == s2->matrix);
}

/*
 * ht_insert_indices()
 *
 * Insert value for (i,j) into the hashtable
 *
 * Parameters:
 *    i,j - indices to build insertion key
 *    matrix - which matrix (I,D,C) to insert for
 *    value - value to insert for the key
 *
 * Return value:
 *    None.
 */
static void ht_insert_indices(unsigned int i, unsigned int j, char matrix, int value)
{
    tuple3_t key;
    unsigned int uvalue;

    key.i = i;
    key.j = j;
    key.matrix = matrix;
    uvalue = value;
    ht_insert(&key, &uvalue);
}

/*
 * ht_lookup_indices()
 *
 * Get the value for (i,j) from the hashtable
 *
 * Parameters:
 *     i,j - indices to build key for lookup
 *     matrix - which matrix (I,D,C) to lookup for
 *     pvalue - (OUTPUT) value for key, only set if TRUE returned
 * 
 * Return value:
 *     TRUE if found, FALSE otherwise
 */
static bool ht_lookup_indices(unsigned int i, unsigned int j, char matrix, int *pvalue)
{
    tuple3_t key;
    unsigned int *pval;

    key.i = i;
    key.j = j;
    key.matrix = matrix;
    pval = (unsigned int *)ht_lookup(&key);
    if (pval)
    {
        *pvalue = *pval;
        return TRUE;
    }
    else
        return FALSE;
}

void modded_three_frame(char* dnaSequence, char* proteinSequence, int C[][strlen(proteinSequence) + 1])
{
    int N = strlen(dnaSequence), M = strlen(proteinSequence);
    int gep = 2, gop = 3, frameshift_penalty = 4;
    int I[N][M + 1], D[N][M + 1];
    int TI[N][M + 1], TD[N][M + 1], TC[N][M + 1];
    unsigned int i,j;

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

    C[0][0] = 0; // this
    TC[0][0] = 0; // this
    for(j = 1; j < M + 1; j++) {
        C[0][j] = 0; // this
        C[j][0] = 0; // and this is already initialized to 0, I think this can be removed (source: initialization step)
        if(!ht_lookup_indices(1, j, 'C', &C[1][j])) {
            C[1][j] = max_of_three(I[1][j], D[1][j], C[0][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 1)));
            ht_insert_indices(1, j, 'C', C[1][j]);
        }
        if(!ht_lookup_indices(2, j, 'C', &C[2][j])) {
            C[2][j] = max_of_two(I[2][j], C[0][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 2)) - frameshift_penalty);
            ht_insert_indices(2, j, 'C', C[2][j]);
        }
        if(!ht_lookup_indices(3, j, 'C', &C[3][j])) {
            C[3][j] = max_of_two(I[3][j], C[1][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 3)) - frameshift_penalty);
            ht_insert_indices(3, j, 'C', C[3][j]);
        }
        if(!ht_lookup_indices(4, j, 'C', &C[4][j])) {
            C[4][j] = max_of_four(I[4][j], D[4][j], C[1][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 4)), C[2][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 4)) - frameshift_penalty);
            ht_insert_indices(4, j, 'C', C[4][j]);
        }

        TC[0][j] = 0; // this
        TC[j][0] = 0; // and this should also be deletable
        if(!ht_lookup_indices(1, j, 'C' + 3, &TC[1][j])) {
            TC[1][j] = C[1][j] == I[1][j] ? -2 : (C[1][j] == D[1][j] ? -1 : (C[1][j] == C[0][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 1)) ? 1 : 0));
            ht_insert_indices(1, j, 'C' + 3, TC[1][j]);
        }
        if(!ht_lookup_indices(2, j, 'C' + 3, &TC[2][j])) {
            TC[2][j] = C[2][j] == I[2][j] ? -2 : (C[2][j] == C[0][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 2)) - frameshift_penalty ? 2 : 0);
            ht_insert_indices(2, j, 'C' + 3, TC[2][j]);
        }
        if(!ht_lookup_indices(3, j, 'C' + 3, &TC[3][j])) {
            TC[3][j] = C[3][j] == I[3][j] ? -2 : (C[3][j] == C[1][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, 3)) - frameshift_penalty ? 2 : 0);
            ht_insert_indices(3, j, 'C' + 3, TC[3][j]);
        }
        if(!ht_lookup_indices(4, j, 'C' + 3, &TC[4][j])) {
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
            ht_insert_indices(4, j, 'C' + 3, TC[4][j]);
        }
    }

    // Change all negative values in C matrix to 0
    for(i = 0; i < N; i++) {
        for(j = 0; j < M + 1; j++) {
            if(C[i][j] < 0)
                C[i][j] = 0;
        }
    }

    // Matrix Filling
    for(i = 0; i < N; i++) {
        for(j = 1; j < M + 1; j++) {
            if(!ht_lookup_indices(i, j, 'I', &I[i][j])) {
                I[i][j] = max_of_two(I[i][j-1] - gep, C[i][j-1] - gop - gep);
                ht_insert_indices(i, j, 'I', I[i][j]);
            }
            if (i < 4)
                continue;
            if(!ht_lookup_indices(i, j, 'D', &D[i][j])) {
                D[i][j] = max_of_two(D[i-3][j] - gep, C[i-3][j] - gop - gep);
                ht_insert_indices(i, j, 'D', D[i][j]);
            }
            if(!ht_lookup_indices(i, j, 'C', &C[i][j])) {
                C[i][j] = max_of_five(
                    I[i][j],
                    D[i][j],
                    C[i-2][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, i)) - frameshift_penalty,
                    C[i-3][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, i)),
                    C[i-4][j-1] + get_score(proteinSequence[j - 1], get_translated_codon(dnaSequence, i)) - frameshift_penalty
                );
                ht_insert_indices(i, j, 'C', C[i][j]);
            }
        }
    }

    // Traceback Matrix Filling
    for(i = 0; i < N; i++) {
        for(j = 0; j < M; j++) {
            if(!ht_lookup_indices(i, j, 'I' + 3, &TI[i][j])) {
                TI[i][j] = I[i][j] == I[i][j-1] - gep ? 0 : (I[i][j] == C[i][j-1] - gop - gep ? 1 : -999);
                ht_insert_indices(i, j, 'I' + 3, TI[i][j]);
            }
            if(!ht_lookup_indices(i, j, 'D' + 3, &TD[i][j])) {
                TD[i][j] = D[i][j] == D[i-3][j] - gep ? 0 : (D[i][j] == C[i-3][j] - gop - gep ? 1 : -999);
                ht_insert_indices(i, j, 'D' + 3, TD[i][j]);
            }
            if(!ht_lookup_indices(i, j, 'C' + 3, &TC[i][j])) {
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
                ht_insert_indices(i, j, 'C' + 3, TC[i][j]);
            }
        }
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
}

int six_frame(char* dnaSequence, char* proteinSequence) {
    int N = strlen(dnaSequence), M = strlen(proteinSequence);
    int C1[N][M + 1], C2[N][M + 1];
    #ifdef PRINTERS
    printf("First Run:\n");
    #endif
    modded_three_frame(dnaSequence, proteinSequence, C1);

    #ifdef PRINTERS
    printf("Reverse Complement:\n");
    #endif
    reverse_complement(dnaSequence);
    modded_three_frame(dnaSequence, proteinSequence, C2);

    int max1 = -999, max2 = -999;
    int val1, val2;
    int i, j;
    for(i = 0; i < strlen(dnaSequence); i++) {
        for(j = 0; j < strlen(proteinSequence) + 1; j++) {
            val1 = C1[i][j];
            val2 = C2[i][j];

            max1 = val1 > max1 ? val1 : max1;
            max2 = val2 > max2 ? val2 : max2;
        }
    }

    return max_of_two(max1, max2);
}

int main() {
    String dnaSequences[] = {"ATGCG", "ATGCGA", "ATGCGATACGCTTGA", "CTTGGTCCGAAT"};
    String proteinSequences[] = {"MCA", "MR", "MRIR", "LGPL"};
    int i;
    double time_taken, start, end;
    for(i = 0; i < 4; i++) {
        printf("DNA Sequence: %s\n", dnaSequences[i]);
        printf("Protein Sequence: %s\n", proteinSequences[i]);
        start = clock();
        printf("Score: %d\n\n", six_frame(dnaSequences[i], proteinSequences[i]));
        end = clock();
        time_taken = ((double)(end - start)*1e3)/CLOCKS_PER_SEC;
        printf("Run %d time taken: %f ms\n\n", i, time_taken);
        // break;
    }
    return 0;
}