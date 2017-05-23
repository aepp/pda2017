#ifndef TASK2_H_INCLUDED
#define TASK2_H_INCLUDED
/* ^^ these are the include guards */

/* Prototypes for the functions */

/**
 * Task 2 of Assignment 2
 *
 * Parameters:
 *
 * argc             number of cli arguments
 * argv             cli arguments
 * sizeOfRandArray  amount of elements in random array
 * commMode         communication mode (non-blocking = 1, persistent = 2)
 */
void task2(int argc, char* argv[], double sizeOfRandArray, double commMode);

#endif