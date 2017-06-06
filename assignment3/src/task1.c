/************************************************************************/
/* Author: Aleksandr Epp <aleksandr.epp@gmail.com>                      */
/* Matriclenumber: 6002853                                              */
/* Assignment : 3                                                       */
/* Task : 1                                                             */
/*                                                                      */
/* Description:                                                         */
/*                                                                      */
/* Matrix multiplication done with Jacobi method                        */
/*                                                                      */
/************************************************************************/

#include "mpi.h"                                // import of the MPI definitions
#include <stdio.h>                              // import of the definitions of the C IO library
#include <string.h>                             // import of the definitions of the string operations
#include <stdlib.h>                             // to dynamically allocate array
#include <unistd.h>                             // standard unix io library definitions and declarations
#include <errno.h>                              // system error numbers
#include <math.h>                               // for log()

#include "task1.h"                              // include own header file
#include "util.h"                               // include assignment utils

#define MAX_RANDOM 100                          // max random

void task1(int argc, char* argv[], double epsilon, char* matrixAFilePath, char* vectorBFilePath, int vectorBSize)
{
    int myRank,                                 // rank of the process
        theirRank,                              // rank of the process sent the message
        DEFAULT_TAG = 1,                        // tag for messages with min
        numProc,                                // number of processes
        i,                                      // iterator variable
        root = 0,                               // root process
        error,                                  // file open error (if any)
        rowsPerProcessCount,                    // rows amount each process gets
        hasConvergedFlag = 0,                   // flag indicating whether desired convergence level reached or not
        iterationsCount = 0;                    // counter for iterations required tii convergence

    double vectorB[vectorBSize],                // vector b
           matrixA[vectorBSize][vectorBSize],   // matrix A
           vectorX[vectorBSize],                // current vector x
           vectorXPrev[vectorBSize];            // previous vector x to calculate the convergence
// -->   double (*myMatrixARows)[vectorBSize],       // rows of matrix A for own calculation
// -->          *myVectorXValues;                    // own calculated vector x values

    MPI_File fh;                                // file handle for mpi file operations

    // initializing of MPI-Interface
    MPI_Init(&argc, &argv);

    MPI_Status status; //capture status of a MPI_Send

    // get your rank
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // get the number of processes running
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);

    // calculate how many rows each process gets
    rowsPerProcessCount = vectorBSize / numProc;

    // allocate space for own rows of matrix A
// -->    myMatrixARows = malloc(sizeof *myMatrixARows * rowsPerProcessCount);

    double myMatrixARows[rowsPerProcessCount][vectorBSize],       // rows of matrix A for own calculation
           myVectorXValues[rowsPerProcessCount];                    // own calculated vector x values
    // test allocation
//    int j;
//    for (i = 0; i < rowsPerProcessCount; i++) {
//        for (j = 0; j < vectorBSize; j++) {
//            myMatrixARows[i][j] = ((double)rand()/(double)RAND_MAX);
//            printf("%lf, ", myMatrixARows[i][j]);
//        }
//    }

    // allocate space for own calculated vector x values
// -->    myVectorXValues = malloc(rowsPerProcessCount);

    // if I'm root
    if (myRank == root){
        // print the parameters
//        printf("Path to matrix A: %s\n", matrixAFilePath);
//        printf("Path to Vector b: %s\n", vectorBFilePath);
//        printf("Epsilon: %lf\n", epsilon);

        // open matrix A file
        error = MPI_File_open(MPI_COMM_SELF, matrixAFilePath, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
        if (error){
            printf("Error while opening matrix A file, make sure the path is correct.\n");
            exit(1);
        };
        // read matrix A from file
        MPI_File_read(fh, matrixA, MPI_File_get_size(fh, 0)/sizeof(double), MPI_DOUBLE, &status);
        /// close matrix A file
        MPI_File_close(&fh);
    }

    // distribute matrix A rows to all processes
    MPI_Scatter(matrixA, vectorBSize * rowsPerProcessCount, MPI_DOUBLE, myMatrixARows, vectorBSize * rowsPerProcessCount, MPI_DOUBLE, root, MPI_COMM_WORLD);

    // test retrieval
//    if(myRank == 2){
//        int j;
//        for (i = 0; i < rowsPerProcessCount; i++) {
//            for (j = 0; j < vectorBSize; j++) {
//                printf("%lf, ", myMatrixARows[i][j]);
//            }
//        }
//        printf("\n");
//    }
    // open vector b file
    error = MPI_File_open(MPI_COMM_SELF, vectorBFilePath, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if (error){
        printf("Error while opening vector b file, make sure the path is correct.\n");
        exit(1);
    };
    // read vector b from file
    MPI_File_read(fh, vectorB, vectorBSize, MPI_DOUBLE, &status);
    /// close vector b file
    MPI_File_close(&fh);

    // test vector b read
//    if(myRank == 2){
//        for (i = 0; i < vectorBSize; i++) {
//            printf("%lf, ", vectorB[i]);
//        }
//        printf("\n");
//    }
    // initialize vector x with 0
    fillWith(vectorX, vectorBSize, 0);
    // initialize vector x with 0
    fillWith(vectorXPrev, vectorBSize, 0);

    // apply Jacobi method until the desired convergence reached
    while (hasConvergedFlag == 0){
        iterationsCount++;

        // calculate new own vector x values
        jacobiIterativeRule(
            vectorBSize,
            myMatrixARows,
            myVectorXValues,
            vectorB,
            vectorXPrev,
            rowsPerProcessCount,
            myRank
        );

        // check new calculated vector x
//        for(i = 0; i < rowsPerProcessCount; i++){
//            printf("myX after: %lf ", myVectorXValues[i]);
//        }
//        printf("\n");

        // collect all new vector x values in vector x
        MPI_Allgather(myVectorXValues, rowsPerProcessCount, MPI_DOUBLE, vectorX, rowsPerProcessCount, MPI_DOUBLE, MPI_COMM_WORLD);

        // set a flag whether desired convergence level reached or not
        hasConvergedFlag = hasConverged(vectorX, vectorXPrev, vectorBSize, epsilon);

        // set previous vector x = current vector x for next iteration
        memcpy(vectorXPrev, vectorX, vectorBSize * sizeof(double));
    }

    // if I'm root, save the results to file
    if (myRank == root){
        // open or create file to save results
        error = MPI_File_open(MPI_COMM_SELF, "result", MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
        // check for errors
        if (error){
            printf("Error while opening file to save results.\n");
            exit(1);
        };
        // write results to file
        error = MPI_File_write(fh, vectorX, vectorBSize, MPI_DOUBLE, &status);
        // check for errors
        if (error){
            printf("Error while writing results to file.\n");
            exit(1);
        };
        // close file
        MPI_File_close(&fh);

        // print the result
        printf("Required iterations: %d \n", iterationsCount);
        printf("The resulting x vector to the equation Ax = b is:\n");
        printf("|-----------------------|\n");
        for (i = 0; i < vectorBSize; i++){
            printf("|\t%lf\t|\n", vectorX[i]);
        }
        printf("|-----------------------|\n");
    }


//    // free myMatrixARows
//    for (i = 0; i < rowsPerProcessCount; i++) {
//      free(myMatrixARows[i]);
//    }
//    free(myMatrixARows);

    // finalizing MPI interface
    MPI_Finalize();
}
