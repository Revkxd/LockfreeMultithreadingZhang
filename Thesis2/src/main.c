#include "lockfreeseq.c"
// #include "recursive_seq.c"
// #define ITER
#ifdef ITER
#include "sequential.c"
#endif

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

    // String dnaSequences[] = {"ATTGACAACCGCGTCCGCCGC","ATTGACAACCGCGTCCGCCGCCGCTTCAAGGGCCAGTACTTGATGCCCAACATTGGCTACGGCTCCAACAAGCGCACCCGCCACATGTTGCCCACCGGCT", "GCTACGTCCGCTCCTCCATGTCCTTGTCCGGCTACATGCCCCCCTTGTGCGACCCCAAGGACGGCCACTTGTTGTTGGACGGCGGCTACGTCAACAACT", "GAGCCCACCTCCGAGATTTTGCAGAACCCCGCCCGCGTCTTGCGCCAGCAGTTGAAGGTCTTGTCCGTCATTGACGGCCAGTCCTACGAGCCCTTGAAGG", "CCCGGCGCCGGCTCCGGCCACGGCCACGGCCCCAACGGCGGCTCCAACTCCTCCTCCTGCACCCCCCCCTCCTCCAACCCCCACATTACCGGCTACGTCG"};
    // String proteinSequences[] = {"IDNRVRR","IDNRVRRRFKGQYLMPNIGYGSNKRTRHMLPTGF", "RYVRSSMSLSGYMPPLCDPKDGHLLLDGGYVNNL", "EPTSEILQNPARVLRQQLKVLSVIDGQSYEPLKD", "PGAGSGHGHGPNGGSNSSSCTPPSSNPHITGYVD"};
    // int correctScores[] = {35, 179, 177, 162, 191};
    // int score;

    // int i;
    // for(i = 0; i < 5; i++) {
        // init_hash_table();
        // printf("DNA Sequence: %s\n", dnaSequences[i]);
        // printf("Protein Sequence: %s\n", proteinSequences[i]);
        // start = clock();
        // printf("Score: %d\n\n", six_frame(dnaSequences[i], proteinSequences[i]));
        // end = clock();
        // time_taken = (double)(end - start)*1e3 / CLOCKS_PER_SEC;
        // printf("Run %d time taken: %f ms\n\n", i, time_taken);
        // break;
    // }

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