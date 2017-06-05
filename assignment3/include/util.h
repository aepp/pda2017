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
void fillWithRandomDouble(double* array, int size, double maxRandom, int rank);

/**
 * fill array with zeroes
 *
 * array        array to fill
 * size         size of the array
 * val          value to fill with
 */
void fillWith(double* array, int size, double val);

/**
 * generate random double value
 *
 * maxRandom    maximum random value
 *
 * return       random double value which is at most [max]
 */
double generateRandomDouble(double max);

/**
 * jacobi iterative rule to calculate vector x values
 *
 * myMatrixARows            own matrix A rows received from root
 * myVectorXValues          array to store new own calculated vector x value
 * vectorB                  vector b
 * vectorXPrev              previous vector x values
 * vectorBSize              size of vector b
 * rowsPerProcessCount      amount of matrix A rows each process gets for calculation
 * myRank                   own rank
 */
void jacobiIterativeRule(
    int vectorBSize,
    double myMatrixARows[][vectorBSize],
    double *myVectorXValues,
    double *vectorB,
    double *vectorXPrev,
    int rowsPerProcessCount,
    int myRank);

/**
 * check whether desired convergence level reached or not
 *
 * vectorX          vector x calculated in current iteration
 * vectorXPrev      vector x calculated in previous iteration
 * vectorBSize      size of vector x
 * epsilon          desired convergence level
 *
 * return           1 in converged to desired epsilon, 0 otherwise
 */
int hasConverged(double *vectorX, double *vectorXPrev, int vectorBSize, double epsilon);

#endif