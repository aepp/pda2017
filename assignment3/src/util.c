#include <stdio.h>                          // import of the definitions of the C IO library
#include <string.h>                         // import of the definitions of the string operations
#include <unistd.h>                         // standard unix io library definitions and declarations
#include <errno.h>                          // system error numbers
#include <stdlib.h>                         // for random()
#include <time.h>                           // to seed random generator
#include <math.h>                           // for math functions
#include "mpi.h"                            // import of the MPI definitions

#include "util.h"                           // include own header file

void fillWithRandomDouble(double* array, int size, double maxRandom, int rank)
{
    srandom((unsigned)time(NULL) * rank);
    int i;
    for(i = 0; i < size; i++) {
        // assign random int to each array position
        array[i] = generateRandomDouble(maxRandom);
    }
}

void fillWith(double* array, int size, double val)
{
    int i;
    for(i = 0; i < size; i++) {
        // assign random int to each array position
        array[i] = val;
    }
}

double generateRandomDouble(double max){
    return ((double)random() / RAND_MAX) * max;
}

void jacobiIterativeRule(
    int vectorBSize,
    double myMatrixARows[][vectorBSize],
    double *myVectorXValues,
    double *vectorB,
    double *vectorXPrev,
    int rowsPerProcessCount,
    int myRank
){
    int row, col,                       // iteration variables
        diagonalElemPos;                // diagonal element position
    double l,                           // lower triangle
           u;                           // upper triangle

    for (row = 0; row < rowsPerProcessCount; row++){
        // determine diagonal element position
        diagonalElemPos = rowsPerProcessCount * myRank + row;
        // reset temporary sums for the next iteration
        l = 0;
        u = 0;

        // left sum in the element-based approach
        for (col = 0; col < diagonalElemPos; col++){
            l += (myMatrixARows[row][col] * vectorXPrev[col]);
        }
        // right sum in the element-based approach
        for (col = diagonalElemPos + 1; col < vectorBSize; col++){
            u += myMatrixARows[row][col] * vectorXPrev[col];
        }
        // calculate new vector x value
        // note: diagonal element always located in myMatrixARows[row][row] in a square matrix
        myVectorXValues[row] = (1.0 / myMatrixARows[row][diagonalElemPos]) * (vectorB[diagonalElemPos] - l - u);
    }
}

int hasConverged(double *vectorX, double *vectorXPrev, int vectorBSize, double epsilon){
    int i;                      // iteration variables
    double distance = 0.0;      // distance between vector x calculated in current and previous iterations

    // iterate through vector x elements
    for (i = 0; i < vectorBSize; i++) {
        // calculate the distance
        distance += fabs(vectorX[i] - vectorXPrev[i]);
    }

    // if desired convergence level reached return 1 (calculation ends)
    if (distance < epsilon){
        return 1;
    }
    // otherwise return 0 (calculation continues)
    return 0;
}