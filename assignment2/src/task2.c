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
#include <string.h>                         // import of the definitions of the string operations
#include <unistd.h>                         // standard unix io library definitions and declarations
#include <errno.h>                          // system error numbers
#include <math.h>                           // for log()

#include "task2.h"                          // include own header file
#include "util.h"                           // include assignment utils

#define ARRAY_SIZE 10                       // size of random array
#define MAX_RANDOM 10                       // max random integer

void task2(int argc, char* argv[])
{
    int myRank,                             // rank of the process
        theirRank,                          // rank of the process sent the message
        DEFAULT_TAG = 1,                    // tag for messages with min
        numProc,                            // number of processes
        i,                                  // iterator variable
        root = 0,                           // root process
        randomInts[ARRAY_SIZE];             // random array
    double timeSeq[2],                      // stores execution time of sequential component
           timePar[2],                      // stores execution time of parallel component
           timeCom[2],                      // stores execution time of all communication components
           myResult,                        // result of own integration
           rResult;                         // received integration results

    // initializing of MPI-Interface
    MPI_Init(&argc, &argv);

    // fill with random integers
    fillWithRandomInt(randomInts, ARRAY_SIZE, MAX_RANDOM, myRank);

    // start time measurement when all initialization is done and communication begins
    timeCom[0] = MPI_Wtime();

    MPI_Status status; //capture status of a MPI_Send

    // get your rank
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // get the number of processes running
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);

    // finalizing MPI interface
    MPI_Finalize();
}
