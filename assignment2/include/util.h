#ifndef UTIL_H_INCLUDED
#define UTIL_H_INCLUDED
/* ^^ these are the include guards */

/* Prototypes for the functions */

/**
 * trapezoidal rule for function 1
 *
 * a        left integration limit
 * b        right integration limit
 * n        amount of sub-intervals
 */
double trapezoidalRuleF1(double a, double b, int n);

/**
 * trapezoidal rule for function 2
 *
 * a        left integration limit
 * b        right integration limit
 * n        amount of sub-intervals
 */
double trapezoidalRuleF2(double a, double b, int n);

/**
 * function to integrate no. 1
 */
double f1(double x);

/**
 * function to integrate no. 2
 */
double f2(double x);

/**
 * Fill array with random integers with the maximum value maxRandom
 *
 * array        array to fill with random integers
 * size         size of array to fill
 * maxRandom    highest random number
 * rank         rank of the current process for random seed
 */
void fillWithRandomInt(int* array, int size, int maxRandom, int rank);

/**
 * generate random integer with the maximum value max
 */
int generateRandomInt(int max);

/**
 * comparator function for quick sort
 */
int cmpFunc (const void *a, const void *b);

/**
 * communication in non-blocking mode, task 2a)
 */
void nonBlockingCommunication();

/**
 * persistent communication, task 2b)
 */
void persistentCommunication();

#endif