/************************************************************************/
/* Author: Aleksandr Epp <aleksandr.epp@gmail.com>                      */
/* Matriclenumber: 6002853                                              */
/* Assignment : 5                                                       */
/* Task : 2                                                             */
/*                                                                      */
/* Description:                                                         */
/*                                                                      */
/* It's a worker main file for the distributed Cannon Algorithm         */
/************************************************************************/

#include "mpi.h"                                // import of the MPI definitions
#include <stdio.h>                              // import of the definitions of the C IO library
#include <string.h>                             // import of the definitions of the string operations
#include <stdlib.h>                             // to dynamically allocate array
#include <unistd.h>                             // standard unix io library definitions and declarations
#include <errno.h>                              // system error numbers
#include <math.h>                               // for log()
#include <sys/select.h>
#include <inttypes.h>                   // for strtoumax() string to int conversion

int main(int argc, char* argv[])
{
    int myRank,                 // rank of the process
        numProc,                // number of processes
        i, j, k, l;

    double *matrixA,            // pointer to use for scatterv()
           *matrixB,            // pointer to use for scatterv()
           *matrixC;            // pointer to use for scatterv()

    MPI_Comm intercomm;
    MPI_Status  status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    MPI_Comm_get_parent(&intercomm);

    // get matrix and sub-matrix dimensions from arguments array
    int subMatrixDimension = strtoumax(argv[1], NULL, 10),
        matrixDimension = strtoumax(argv[2], NULL, 10);

    double mySubMatrixA[subMatrixDimension][subMatrixDimension],
           mySubMatrixB[subMatrixDimension][subMatrixDimension],
           mySubMatrixC[subMatrixDimension][subMatrixDimension];

    // define datatype for a single sub-matrix
    MPI_Datatype tmpDatatype, subMatrixDatatype;
    int matrixSizes[2] = {matrixDimension, matrixDimension},
        subMatrixSizes[2] = {subMatrixDimension, subMatrixDimension},
        starts[2] = {0, 0};

    MPI_Type_create_subarray(2, matrixSizes, subMatrixSizes, starts, MPI_ORDER_C, MPI_DOUBLE, &tmpDatatype);
    MPI_Type_create_resized(tmpDatatype, 0, subMatrixDimension * sizeof(double), &subMatrixDatatype);
    MPI_Type_commit(&subMatrixDatatype);

    int sendCountsPerProcess[numProc],          // amount of units per process
        displacementsPerProcess[numProc],       // displacements per process
        tmpDisplacement = 0;

    for(i = 0; i < numProc; i++){
        sendCountsPerProcess[i] = 1;
    }

    for(i = 0; i < numProc; i++){
        displacementsPerProcess[i] = 1;
    }

    for(i = 0; i < sqrt(numProc); i++){
        for(j = 0; j < sqrt(numProc); j++){
            displacementsPerProcess[i * (int)sqrt(numProc) + j] = tmpDisplacement;
//            printf("worker %d\n ", i * (int)sqrt(numProc) + j);
            tmpDisplacement++;
        }
        tmpDisplacement += (subMatrixDimension - 1) * (int)sqrt(numProc);
    }

//    printf("Displacements for %d subMatrixDimension: ", subMatrixDimension);
//    for(i = 0; i < 2; i++){
//        printf("%d, ", displacementsPerProcess[i]);
//    }
//    printf("\n");

    // scatter matrix A
    MPI_Scatterv(
        matrixA,
        sendCountsPerProcess,
        displacementsPerProcess,
        subMatrixDatatype,
        mySubMatrixA,
        subMatrixDimension * subMatrixDimension,
        MPI_DOUBLE,
        //MPI_ROOT,
        0,
        intercomm
    );

    // scatter matrix B
    MPI_Scatterv(
        matrixB,
        sendCountsPerProcess,
        displacementsPerProcess,
        subMatrixDatatype,
        mySubMatrixB,
        subMatrixDimension * subMatrixDimension,
        MPI_DOUBLE,
        //MPI_ROOT,
        0,
        intercomm
    );


    // create communicator topologies
    int dimensionCount = 2,
        dimensionSize[2] = {sqrt(numProc), sqrt(numProc)},
        wrap[2] = {1, 1},
        reorder = 1;

    MPI_Comm matrixACommunicator,
             rowsOfMatrixACommunicator,
             matrixBCommunicator,
             colsOfMatrixBCommunicator;

    // for matrix A
    MPI_Cart_create(MPI_COMM_WORLD, dimensionCount, dimensionSize, wrap, reorder, &matrixACommunicator);

    // for matrix B
    MPI_Cart_create(MPI_COMM_WORLD, dimensionCount, dimensionSize, wrap, reorder, &matrixBCommunicator);

    // create sub-grid for rows of matrix A
    int remainDims[2] = {0, 1}; // keep columns, drop rows

    MPI_Cart_sub(matrixACommunicator, remainDims, &rowsOfMatrixACommunicator);
    // get process ranks in row grid
//    MPI_Comm_rank(rowsOfMatrixACommunicator, &row_rankA);

    // create sub-grid for columns of matrix B
    remainDims[0] = 1; // keep columns
    remainDims[1] = 0; // drop rows

    MPI_Cart_sub(matrixBCommunicator, remainDims, &colsOfMatrixBCommunicator);
    // get process ranks in column grid
//    MPI_Comm_rank(colsOfMatrixBCommunicator, &col_rankB);

    int myRowId = (int)(myRank/sqrt(numProc)),  // get id of my row grid in matrix A
        myColId = myRank%(int)sqrt(numProc),    // get id of my column grid in matrix B
        srcRank,                                // process id
        destRank;

    // printf("rank: %d with row index %d\n", myRank, myRowId);
    // printf("rank: %d with column index %d\n", myRank, myColId);

    // calculate source and destination ranks
    // for initial shift of matrix A
    MPI_Cart_shift(rowsOfMatrixACommunicator, 0, -myRowId, &srcRank, &destRank);
    // and exchange values
    MPI_Sendrecv_replace(
        mySubMatrixA,
        subMatrixDimension * subMatrixDimension,
        MPI_DOUBLE,
        destRank,
        0,
        srcRank,
        0,
        rowsOfMatrixACommunicator,
        &status
    );

    // calculate source and destination ranks
    // for initial shift of matrix B
    MPI_Cart_shift(colsOfMatrixBCommunicator, 0, -myColId, &srcRank, &destRank);

    // and exchange values
    MPI_Sendrecv_replace(
        mySubMatrixB,
        subMatrixDimension * subMatrixDimension,
        MPI_DOUBLE,
        destRank,
        0,
        srcRank,
        0,
        colsOfMatrixBCommunicator,
        &status
    );

    // wait until every process exchanged values
    MPI_Barrier(MPI_COMM_WORLD);

//    printf("A:\n");
//    for (k = 0; k < subMatrixDimension; k++){
//        for (l = 0; l < subMatrixDimension; l++){
//            printf("%.1f, ", mySubMatrixA[k][l]);
//        }
//        printf("\n");
//    }
//    printf("B:\n");
//    for (k = 0; k < subMatrixDimension; k++){
//        for (l = 0; l < subMatrixDimension; l++){
//            printf("%.1f, ", mySubMatrixB[k][l]);
//        }
//        printf("\n");
//    }
//    printf("C (init):\n");
//    for (k = 0; k < subMatrixDimension; k++){
//        for (l = 0; l < subMatrixDimension; l++){
//            printf("%.1f, ", mySubMatrixC[k][l]);
//        }
////        if(mySubMatrixC[j][k] > 999999){
////            printf("j: %d; k: %d\n", j, k);
////        }
//        printf("\n");
//    }
    for (j = 0; j < subMatrixDimension; j++){
        for (k = 0; k < subMatrixDimension; k++){
            mySubMatrixC[j][k] = 0.0; //initialize empty result matrix with zeros
        }
    }
    double tmpSum;
    for(i = 0; i < sqrt(numProc); i++){
        // do multiplication
        for (j = 0; j < subMatrixDimension; j++){
            for (k = 0; k < subMatrixDimension; k++){
                tmpSum = 0.0;
                for (l = 0; l < subMatrixDimension; l++){
                    tmpSum += mySubMatrixA[j][l] * mySubMatrixB[l][k];
                }
                mySubMatrixC[j][k] += tmpSum;
            }
//            printf("\n");
        }
//        if(myRank == 0){
//            for (k = 0; k < subMatrixDimension; k++){
//                for (l = 0; l < subMatrixDimension; l++){
//                    printf("%.1f, ", mySubMatrixC[k][l]);
//                }
//                printf("\n");
//            }
//        }
        MPI_Cart_shift(rowsOfMatrixACommunicator, 0, -1, &srcRank, &destRank);
        MPI_Sendrecv_replace(
            mySubMatrixA,
            pow(subMatrixDimension, 2),
            MPI_DOUBLE,
            destRank,
            0,
            srcRank,
            0,
            rowsOfMatrixACommunicator,
            &status
        );

        MPI_Cart_shift(colsOfMatrixBCommunicator, 0, -1, &srcRank, &destRank);
        MPI_Sendrecv_replace(
            mySubMatrixB,
            pow(subMatrixDimension, 2),
            MPI_DOUBLE,
            destRank,
            0,
            srcRank,
            0,
            colsOfMatrixBCommunicator,
            &status
        );
    }

    // wait until everyone has finished calculations
    MPI_Barrier(MPI_COMM_WORLD);

    // process 0 sends info to master process that calculations are finished
    if (myRank == 0){
        int done = 1;    
        MPI_Request doneRequest;

        MPI_Isend(&done, 1, MPI_INT, 0, 0, intercomm, &doneRequest);
        MPI_Wait(&doneRequest, &status);
    }


//    printf("\n\n");
//    for (i = 0; i < 16; i++){
//        for (j = 0; j < 16; j++){
//            printf("%.1lf ", mySubMatrixC[i][j]);
//            if(j == 15) printf("\n");
//        }
//    }
//    printf("\n\n");

//    if(myRank == 3){
//        sleep(3);
//    }
//    if(myRank == 1){
//        printf("\n\n");
//        for (i = 0; i < 8; i++){
//            for (j = 0; j < 8; j++){
//                printf("%.1lf ", mySubMatrixC[i][j]);
//                if(j == 7) printf("\n");
//            }
//        }
//        printf("\n\n");
//    }
//    if(myRank == 3){
//        printf("\n\n");
//        for (i = 0; i < 8; i++){
//            for (j = 0; j < 8; j++){
//                printf("%.1lf ", mySubMatrixC[i][j]);
//                if(j == 7) printf("\n");
//            }
//        }
//        printf("\n\n");
//    }

    // process 0 collects the results
    MPI_Gatherv(
        mySubMatrixC,
        subMatrixDimension * subMatrixDimension,
        MPI_DOUBLE,
        matrixC,
        sendCountsPerProcess,
        displacementsPerProcess,
        subMatrixDatatype,
        0,
        intercomm
    );

    MPI_Type_free(&subMatrixDatatype);
    MPI_Finalize();

    return 0;
}
