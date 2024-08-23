# Integrating Lock-Free Multithreading to Modified Zhang’s Six-Frame Alignment Algorithm
This is an exploratory study based on the Six-Frame modification of Zhang's Frame Alignment algorithm aimed at comparing the speed and accuracy of an implementation using multithreading (using the pthread library) with lock-free hash tables. Three versions of the algorithm were implemented, where a bottom-up version was implemented as a baseline for accuracy and 2 top-down approaches (single-threaded, multithreaded). The dataset used was the sequence of Drosophila melanogaster (fruit fly) retrieved from the UniProt website. In testing, 4 unique DNA sequences were compared against a set of 1000 protein sequences each with the score and runtime of the different versions being compared against each other.

This is an exploratory study aimed at exploring a new implementation approach for Zhang’s six-frame alignment. The group wanted to investigate its feasibility, especially considering the potential challenges with multiple matrices. Contrary to the theoretical advantages of multithreading, the lock-free version using pthreads and increasing thread counts was consistently slower compared to the recursive version. This slower performance can be attributed to the additional computation at the start which is based on a random index and the atomic operation used in ht_insert. Future work may explore optimizing the implementation of multithreaded approach or identifying scenarios where its benefits can be fully realized.

# Adviser
### Mr. Roger Luis Uy (roger.uy@dlsu.edu.ph)

# Members
- Kalaw, Stacy Selena (stacy_selena_kalaw@dlsu.edu.ph)
- Meneses, Alyssa (alyssa_meneses@dlsu.edu.ph)
- Sy, Richard (richard_t_sy@dlsu.edu.ph)
- Yu, Cedric Leopold (cedric_leopold_yu@dlsu.edu.ph)

# Relevant Files Information (under Thesis2/src folder)
## Main Files (dynamic programming algorithm codes)
- lockfreeseq.c -> Implementation of Multithreading (using pthread library) with lock-free hash table top-down version
- recursive_seq.c -> Implementation of single-threaded top-down version
- sequential.c -> Latest, fixed version of `seq.c` which is the bottom-up implementation
- seq.c -> Old, wrong implementation, ignore this file.

## Helper Files
- constants.h -> constants definitions and some extra user defined data types are stored here
- tables.h -> Stores the codon translation table and blosum62 scoring matrix
- maxfunctions.c -> function implementations of `max()` function to allow for different argument counts
- converters.c -> main helper file storing DNA to Protein conversion helpers, and retrieves score from the blosum matrix
- splitlfht.c -> Implementation of lock-free hash table used in `lockfreeseq.c` and `recursive_seq.c`

# Steps to Run
1. Clone the repository in a Linux system (we use the pthread.h library which is exclusive to Linux)
2. cd into the Thesis2/src folder as the main code used in the paper is contained in that folder
3. Since there are 3 versions of code, you compile whichever version you want from the **Main Files** section
4. The next steps require you to be in your command line in `Thesis2/src` folder
5. To compile the sequential version, run `gcc remake.c`
6. To compile the top-down versions, run `gcc -O3 <filename of top-down version>` (ex. lockfreeseq.c)
7. You may now execute the program with `./a.out <DNA sequence> <Protein sequence>` which some examples are provided in the `testcases.txt` file

# Repository Sources
- [Lockfree Hash Table](https://github.com/stivalaa/paralleldp)

# Side Notes
- Makefile is not configured properly since the usage was dropped early into the development
- The main files are crammed with all the functions to have proper debug outputs which could've been avoided if a Makefile was configured
- Most of the execution for tests were done on Google Colab which is why the Makefile was dropped to have easier code migration to Colab after finishing coding