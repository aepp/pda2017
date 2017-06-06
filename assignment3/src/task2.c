/************************************************************************/
/* Author: Aleksandr Epp <aleksandr.epp@gmail.com>                      */
/* Matriclenumber: 6002853                                              */
/* Assignment : 3                                                       */
/* Task : 2                                                             */
/*                                                                      */
/* Description:                                                         */
/*                                                                      */
/* Matrix multiplication done with Jacobi method.                       */
/* Same Jacobi method as implemented in task 1 but tis time processes   */
/* read columns and not rows from the input file(s).                    */
/*                                                                      */
/************************************************************************/

#include "mpi.h"                                // import of the MPI definitions
#include <stdio.h>                              // import of the definitions of the C IO library
#include <string.h>                             // import of the definitions of the string operations
#include <stdlib.h>                             // to dynamically allocate array
#include <unistd.h>                             // standard unix io library definitions and declarations
#include <errno.h>                              // system error numbers
#include <math.h>                               // for log()

#include "task2.h"                              // include own header file
#include "util.h"                               // include assignment utils

#define MAX_RANDOM 100                          // max random

void task2(int argc, char* argv[], double epsilon, char* matrixAFilePath, char* vectorBFilePath, int vectorBSize)
{
    int myRank,                                 // rank of the process
        theirRank,                              // rank of the process sent the message
        DEFAULT_TAG = 1,                        // tag for messages with min
        numProc,                                // number of processes
        i,                                      // iterator variable
        root = 0,                               // root process
        error,                                  // file open error (if any)
        colsPerProcessCount,                    // rows amount each process gets
        hasConvergedFlag = 0,                   // flag indicating whether desired convergence level reached or not
        iterationsCount = 0,                    // counter for iterations required tii convergence
        displacement;                           // position for each process to read the file

    double matrixA[vectorBSize][vectorBSize],   // matrix A
           vectorX[vectorBSize],                // current vector x
           vectorXPrev[vectorBSize];            // previous vector x to calculate the convergence
// -->    double *myVectorXValues;

    MPI_File fh;                                // file handle for mpi file operations

    MPI_Datatype stripes;                       // new datatype to read columns from a matrix

    // initializing of MPI-Interface
    MPI_Init(&argc, &argv);

    MPI_Status status; //capture status of a MPI_Send

    // get your rank
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // get the number of processes running
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);

    // calculate how many rows each process gets
    colsPerProcessCount = vectorBSize / numProc;

    // set displacement
    displacement = myRank * sizeof(double) * colsPerProcessCount;

    // define new data type
    MPI_Type_vector(vectorBSize, colsPerProcessCount, vectorBSize, MPI_DOUBLE, &stripes);
    // and commit it
    MPI_Type_commit(&stripes);

    // allocate space for own columns of matrix A, vector B and vector X
    double myMatrixAColumns[vectorBSize][colsPerProcessCount],  // own matrix A columns
           myVectorB[colsPerProcessCount],                      // own part of vector b
           myVectorXValues[colsPerProcessCount];                // own calculated vector x values

// -->    myVectorXValues = malloc(colsPerProcessCount);

    // open matrix A file
    error = MPI_File_open(MPI_COMM_WORLD, matrixAFilePath, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if (error){
        printf("Error while opening matrix A file, make sure the path is correct.\n");
        exit(1);
    };
    // set file view to be able to read stripes from the matrix
    error = MPI_File_set_view(fh, displacement, MPI_DOUBLE, stripes, "native", MPI_INFO_NULL);
    if (error){
        printf("Error while setting file view.\n");
        exit(1);
    };
    // read stripes from matrix A file (split collective, shared file pointer, doesn't work)
//    error = MPI_File_read_ordered_begin(fh, myMatrixAColumns, colsPerProcessCount * vectorBSize, MPI_DOUBLE);
//    error = MPI_File_read_ordered_end(fh, &myMatrixAColumns, &status);

    // read stripes from matrix A file (collective, individual file pointer, works)
    error = MPI_File_read_all(fh, myMatrixAColumns, colsPerProcessCount * vectorBSize, MPI_DOUBLE, &status);

    // close matrix A file
    MPI_File_close(&fh);

    // test split collective read from file with individual pointers
//    if(myRank == 1){
//        int j;
//        for(i = 0; i < vectorBSize; i++){
//            for(j = 0; j < colsPerProcessCount; j++){
//                printf("%d: %lf\n  ", myRank, myMatrixAColumns[i][j]);
//            }
//        }
//        printf("\n");
//    }

    // open vector b file
    error = MPI_File_open(MPI_COMM_WORLD, vectorBFilePath, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if (error){
        printf("Error while opening vector b file, make sure the path is correct.\n");
        exit(1);
    };
    // read from vector b file (split collective, shared file pointer, works)
    MPI_File_read_ordered(fh, myVectorB, colsPerProcessCount, MPI_DOUBLE, &status);
    // close vector b file
    MPI_File_close(&fh);

    // test read from vector b file
//    if(myRank == 3){
//        for(i = 0; i < colsPerProcessCount; i++){
//            printf("%d: %lf,  ", myRank, myVectorB[i]);
//        }
//        printf("\n");
//    }

    // initialize vector x with 0
    fillWith(vectorX, vectorBSize, 20);
    // initialize vector x with 0
    fillWith(vectorXPrev, vectorBSize, 0);

    double tmpSumLocal, sumGlobal;  // partial and complete sums of iterative rule
    int reducerRank,                // process rank which receives all parts of x through reduce operation
        row, col,                   // iterator variables
        diagElemProcRank,           // process containing diagonal element for the current row
        diagElemPos;                // position of diagonal element in the current row

    // apply Jacobi method until the desired convergence reached
    while (hasConvergedFlag == 0){
        iterationsCount++;

        // calculate new own vector x values
        for (row = 0; row < vectorBSize; row++){
            reducerRank = row / colsPerProcessCount;
            // reset temp sum
            tmpSumLocal = 0;
            // go through own matrix A columns
            for (col = 0; col < colsPerProcessCount; col++){
                // determine diagonal element position
                diagElemPos = myRank * colsPerProcessCount + col;

                // if current position is not the diagonal element
                if (row != diagElemPos){
                    // update temp sum via iterative rule
                    tmpSumLocal += myMatrixAColumns[row][col] * vectorXPrev[diagElemPos];
                } else {
                    // otherwise remember process holding the diagonal element
                    diagElemProcRank = col;
                }
            }

            MPI_Reduce(&tmpSumLocal, &sumGlobal, 1, MPI_DOUBLE, MPI_SUM, reducerRank, MPI_COMM_WORLD);

            if(myRank == reducerRank){
                myVectorXValues[diagElemProcRank] =
                    (1 / myMatrixAColumns[row][diagElemProcRank]) * (myVectorB[diagElemProcRank] - sumGlobal);
            }
        }

        // check new calculated vector x
//        for(i = 0; i < rowsPerProcessCount; i++){
//            printf("myX after: %lf ", myVectorXValues[i]);
//        }
//        printf("\n");

        // collect all new vector x values in vector x
        MPI_Allgather(myVectorXValues, colsPerProcessCount, MPI_DOUBLE, vectorX, colsPerProcessCount, MPI_DOUBLE, MPI_COMM_WORLD);

        // set a flag whether desired convergence level reached or not
        hasConvergedFlag = hasConverged(vectorX, vectorXPrev, vectorBSize, epsilon);

        // set previous vector x = current vector x for next iteration
        memcpy(vectorXPrev, vectorX, vectorBSize * sizeof(double));
    }

    if (myRank == root){
        // open or create file to save results
        error = MPI_File_open(MPI_COMM_SELF, "result_task_2", MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
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

    // finalizing MPI interface
    MPI_Finalize();
}
