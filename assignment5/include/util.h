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

#endif