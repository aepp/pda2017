#ifndef TASK1_H_INCLUDED
#define TASK1_H_INCLUDED
/* ^^ these are the include guards */

/* Prototypes for the functions */

/**
 * Task 1 of Assignment 3
 *
 * Parameters:
 *
 * argc                     number of cli arguments
 * argv                     cli arguments
 * epsilon                  epsilon value for Jacobi method
 * matrixAFilePath          file name where matrix is stored
 * vectorBFilePath          file name where vector b is stored
 * vectorBSize              size of vector b
 */
void task1(int argc, char* argv[], double epsilon, char* matrixAFilePath, char* vectorBFilePath, int vectorBSize);

#endif