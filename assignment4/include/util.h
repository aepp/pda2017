#ifndef UTIL_H_INCLUDED
#define UTIL_H_INCLUDED
/* ^^ these are the include guards */

/* Prototypes for the functions */

/**
 * fill array with random double values
 *
 * array        array to fill
 * size         size of the array
 * maxRandom    maximum random value
 * rank         own process rank to make values real randoms
 */
void fillWithRandomInt(int* array, int size, int maxRandom, int rank);

/**
 * generate random int value
 *
 * maxRandom    maximum random value
 *
 * return       random int value which is at most [max]
 */
int generateRandomInt(int max);

/**
 * compare function for ascending sort
 */
int cmpFuncASC (const void *a, const void *b);

/**
 * compare function for descending sort
 */
int cmpFuncDESC (const void *a, const void *b);

/**
 * odd-even sort from task 2
 *
 * additional parameter sortOrder for even/odd rows;
 * get rank and comm size from communicator
 */
void oddEvenTranspositionSort(int *myRandomInts, MPI_Comm comm, int sizeOfRandArray, int sortOrder);

/**
 * used internally by odd-even sort from task 2
 */
void getMyNewArray(int *myRandomInts, int *theirRandomInts, int sizeOfRandArray, int partner, int myRank, int sortOrder);

#endif