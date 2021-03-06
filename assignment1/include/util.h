#ifndef UTIL_H_INCLUDED
#define UTIL_H_INCLUDED
/* ^^ these are the include guards */

/* Prototypes for the functions */

/**
 * Fill array with random integers with the maximum value maxRandom
 *
 * array        array to fill with random integers
 * size         size of array to fill
 * rank         rank of the current process for random seed
 */
void fillWithRandomInt(int* array, int size, int maxRandom, int rank);

// get minimum value of array
int getMin(int* array, int size);

// get maximum value of array
int getMax(int* array, int size);

// get sum of array elements
int getSum(int* array, int size);

// generate random integer with the maximum value max
int generateRandomInt(int max);

#endif