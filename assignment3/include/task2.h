#ifndef TASK2_H_INCLUDED
#define TASK2_H_INCLUDED
/* ^^ these are the include guards */

/* Prototypes for the functions */

/**
 * Task 2 of Assignment 3
 *
 * Parameters:
 *
 * argc                     number of cli arguments
 * argv                     cli arguments
 * epsilon                  epsilon value for Jacobi method
 * matrixFileName           file name where matrix is stored
 * vectorBFileName          file name where vector b is stored
 * vectorBSize              size of vector b
 */
void task2(int argc, char* argv[], double epsilon, char* matrixFileName, char* vectorBFileName, int vectorBSize);

#endif