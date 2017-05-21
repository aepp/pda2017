#ifndef UTIL_H_INCLUDED
#define UTIL_H_INCLUDED
/* ^^ these are the include guards */

/* Prototypes for the functions */

/**
 * trapezoidal rule for function 1
 */
double trapezoidalRuleF1(double a, double b, double s, int n);

/**
 * trapezoidal rule for function 2
 */
double trapezoidalRuleF2(double a, double b, double s, int n);

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

#endif