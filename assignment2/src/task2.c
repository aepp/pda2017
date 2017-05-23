/************************************************************************/
/* Author: Aleksandr Epp <aleksandr.epp@gmail.com>                      */
/* Matriclenumber: 6002853                                              */
/* Assignment : 2                                                       */
/* Task : 2                                                             */
/*                                                                      */
/* Description:                                                         */
/*                                                                      */
/************************************************************************/

#include "mpi.h"                            // import of the MPI definitions
#include <stdio.h>                          // import of the definitions of the C IO library
#include <stdlib.h>                         // to use exit() function
#include <string.h>                         // import of the definitions of the string operations
#include <unistd.h>                         // standard unix io library definitions and declarations
#include <errno.h>                          // system error numbers
#include <math.h>                           // for log()

#include "task2.h"                          // include own header file
#include "util.h"                           // include assignment utils

#define MAX_RANDOM 10                       // max random integer

void task2(int argc, char* argv[], double sizeOfRandArray, double commMode)
{
    int myRank,                             // rank of the process
        theirRank,                          // rank of the process sent the message
        DEFAULT_TAG = 1,                    // tag for messages with min
        numProc,                            // number of processes
        i,                                  // iterator variable
        root = 0,                           // root process
        randomInts[(int)sizeOfRandArray],        // random array
        *finalResult;                        // final result of odd-even sort

    // initializing of MPI-Interface
    MPI_Init(&argc, &argv);

    MPI_Status status; //capture status of a MPI_Send

    // get your rank
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // get the number of processes running
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);

    // fill with random integers
    fillWithRandomInt(randomInts, (int)sizeOfRandArray, MAX_RANDOM, myRank);

    // sort own array of random integers
    qsort(randomInts, (int)sizeOfRandArray, sizeof(int), cmpFunc);

    // decide how to communicate
    switch((int)commMode){
        case 1: // non-blocking
            nonBlockingCommunication(randomInts, myRank, numProc, (int)sizeOfRandArray);
            break;
        case 2: // persistent
            persistentCommunication(randomInts, myRank, numProc, (int)sizeOfRandArray);
            break;
        default: // unknown communication mode
            exit(1);
    }

    // root process allocates memory for sorted result
    if(myRank == root){
        finalResult = malloc(numProc * (int)sizeOfRandArray * sizeof(int));
    }
    // root process  collect results from all processes
    MPI_Gather(&randomInts, (int)sizeOfRandArray, MPI_INT, finalResult, (int)sizeOfRandArray, MPI_INT, root, MPI_COMM_WORLD);

    // root process prints the sorted result
    if(myRank == root){
        printf("\nSorted result:\n");
        for(i = 0; i < numProc * (int)sizeOfRandArray; i++){
            printf("%d, ", finalResult[i]);
            if ((i + 1) % (int)sizeOfRandArray == 0){
                printf("\n");
            }
        }
    }
    // finalizing MPI interface
    MPI_Finalize();
}
