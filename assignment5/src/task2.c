/************************************************************************/
/* Author: Aleksandr Epp <aleksandr.epp@gmail.com>                      */
/* Matriclenumber: 6002853                                              */
/* Assignment : 5                                                       */
/* Task : 2                                                             */
/*                                                                      */
/* Description:                                                         */
/*                                                                      */
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

void task2(int argc, char* argv[])
{
    int myRank,                                 // rank of the process
        theirRank,                              // rank of the process sent the message
        DEFAULT_TAG = 1,                        // tag for messages with min
        numProc,                                // number of processes
        i,                                      // iterator variable
        root = 0,                               // root process
        error,                                  // file open error (if any)
        rowsPerProcessCount;                    // amount of image rows per process



    MPI_File fh;                                // file handle for mpi file operations

    // initializing of MPI-Interface
    MPI_Init(&argc, &argv);

    MPI_Status status; //capture status of a MPI_Send

    // get your rank
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // get the number of processes running
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);

    // finalizing MPI interface
    MPI_Finalize();
}
