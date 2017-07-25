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
#include <sys/select.h>

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
        rowsPerProcessCount,                    // amount of image rows per process

        fd = 0,                                 // file descriptor
        selectResult,                           // select result
        readResult,                             // result of the read() function
        maxMultLimit = 5,                       // maximum amount of parallel multiplications
        activeMultCount = 0,                    // amount of currently running multiplications
        workerCounts[maxMultLimit],             // amount of workers per multiplication
        matricesDimensions[maxMultLimit],       // dimensions of original matrices per multiplication
        quit = 0;                               // flag indicating whether user entered q or not


    fd_set readfds;                             // array of sockets to be checked for readability
    struct timeval timeout;                     // timeout for select()
    char buf[11];                               // buffer for input from stdin
    MPI_File fh;                                // file handle for mpi file operations
    MPI_Comm intercomms[maxMultLimit];          // array of intercommunicators
    MPI_Request openRequests[maxMultLimit];     // array of requests
    MPI_Datatype blockTypes[maxMultLimit];      // array of sub matrix datatypes

    // initializing of MPI-Interface
    MPI_Init(&argc, &argv);

    MPI_Status status; //capture status of a MPI_Send

    // get your rank
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // get the number of processes running
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    // endless until user types quit
    while (quit == 0){
        FD_ZERO(&readfds);
        FD_SET(fd, &readfds);

        // set timeout for user input to 5 seconds
        timeout.tv_sec = 5;     // seconds
        timeout.tv_usec = 0;    // microseconds

        // ask user to input the action
        printf("Please enter (s)tart or (q)uit:\n");
        // start select
        selectResult = select(1, &readfds, NULL, NULL, &timeout);

        // if input not empty
        if(selectResult != 0){
            // read input
            memset((void *) buf, 0, 11);
            readResult = read(fd, (void *) buf, 10);

            // if user entered s (=start)
            if (buf[0] == 's'){
                if (activeMultCount == maxMultLimit){
                    printf("You have reached the maximum number of parallel multiplications. Wait for one to finish.\n");
                } else {
                    // increase the counter of currently running multiplications
                    activeMultCount++;
                    // start actual multiplication
                    startMultiplications(
                        intercomms,
                        openRequests,
                        blockTypes,
                        workerCounts,
                        activeMultCount,
                        matricesDimensions
                    );
                }
            } else if(buf[0] == 'q'){ // if user entered q (=quit)
                quit = 1;
            }
        }
        else { // if no input within input timeout
            // and any multiplication still running
//            printf("activeMultCount: %d\n", activeMultCount);
            if (activeMultCount > 0){
                printf("%d multiplication(s) still running...\n", activeMultCount);
                // check how many calculations are still running
                checkForRunningMultiplications(
                    intercomms,
                    openRequests,
                    blockTypes,
                    &activeMultCount,
                    matricesDimensions,
                    workerCounts
                );
//                activeMultCount--;
//                quit = 1;
            }
            else { // otherwise print quit message
                printf("No running multiplications...\n");
            }
        }
    }

    waitForRunningMultiplications(
        intercomms,
        openRequests,
        blockTypes,
        &activeMultCount,
        matricesDimensions,
        workerCounts
    );

    // finalizing MPI interface
    MPI_Finalize();
}
