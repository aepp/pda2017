/************************************************************************/
/* Author: Aleksandr Epp <aleksandr.epp@gmail.com>                      */
/* Matriclenumber: 6002853                                              */
/* Assignment : 2                                                       */
/* Task : 1                                                             */
/*                                                                      */
/* Description:                                                         */
/*                                                                      */
/* User defines integration limits and each process integrate a         */
/* subinterval in those limits. Root process broadcasts user input to   */
/* other processes. After that each process calculates own integration  */
/* limits and integrate the desired function. Integration results will  */
/* be sent form each process to partner process using butterfly         */
/* communication algorithm                                              */
/*                                                                      */
/************************************************************************/

#include "mpi.h"                            // import of the MPI definitions
#include <stdio.h>                          // import of the definitions of the C IO library
#include <string.h>                         // import of the definitions of the string operations
#include <stdlib.h>                         // to dynamically allocate array
#include <unistd.h>                         // standard unix io library definitions and declarations
#include <errno.h>                          // system error numbers
#include <math.h>                           // for log()

#include "task1.h"                          // include own header file
#include "util.h"                           // include assignment utils

//#define MPI_WTIME_IS_GLOBAL 1

void task1(int argc, char* argv[], double a, double b, double n, double funcNumber)
{
    int myRank,                             // rank of the process
        theirRank,                          // rank of the process sent the message
        DEFAULT_TAG = 1,                    // tag for messages with min
        numProc,                            // number of processes
        i,                                  // iterator variable
        root = 0,                           // root process
        numberOfParameters = 4,             // number of params to broadcast
        shift,                              // value to shift for partner rank determination
        bufferSize;                         // size of communication buffer

    double timeSeq[2],                      // stores execution time of sequential component
           timePar[2],                      // stores execution time of parallel component
           timeCom[2],                      // stores execution time of all communication components
           myResult,                        // result of own integration
           rResult,                         // received integration results
           completeInterval,                // the whole integration interval
           subIntervalLength,               // individual integration interval for each process
           myA,                             // individual left integration limit
           myB,                             // individual right integration limit
           parameters[numberOfParameters],  // store user input for broadcasting here
           rParameters[numberOfParameters], // received parameters
           *buffer;                          // communication buffer

    // store user input
    parameters[0] = a;
    parameters[1] = b;
    parameters[2] = n;
    parameters[3] = funcNumber;

    // initializing of MPI-Interface
    MPI_Init(&argc, &argv);

    // start time measurement when all initialization is done and communication begins
    timeCom[0] = MPI_Wtime();

    MPI_Status status; //capture status of a MPI_Send

    // get your rank
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // get the number of processes running
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);

    // determine buffer size
    bufferSize = sizeof(double) * (numProc + MPI_BSEND_OVERHEAD);
    // allocate memory for communication buffer
    buffer = malloc(bufferSize);

    // and attach it
    MPI_Buffer_attach(buffer, bufferSize);

    // if I'm the first process:
    // tell everyone what to do
    // otherwise:
    // receive info about what to do
    MPI_Bcast(parameters, numberOfParameters, MPI_DOUBLE, myRank, MPI_COMM_WORLD);

    // wait until everyone received their value
    MPI_Barrier(MPI_COMM_WORLD);

    // read received parameters into variables
    a = parameters[0];
    b = parameters[1];
    n = parameters[2];
    funcNumber = parameters[3];

    // determine the complete interval
    completeInterval = b - a;
    // determine single sub-interval size
    subIntervalLength = completeInterval / numProc;
    // determine own left integration limit
    myA = a + subIntervalLength * myRank;
    // determine own right integration limit
    myB = a + subIntervalLength * (myRank + 1);

    // integrate
    // use f1 or f2 depending on user input
    if(funcNumber == 1) {
        myResult = trapezoidalRuleF1(myA, myB, n);
    } else {
        myResult = trapezoidalRuleF2(myA, myB, n);
    }

    // communicate like a beautiful butterfly
    shift = numProc;
    for (i = 0; i < log(numProc)/log(2); i++){
        shift = shift >> 1; // shift for partner rank determination

        // printf("I am process %d and I am sending to process %d\n", myRank, myRank^shift);

        MPI_Bsend(&myResult, 1, MPI_DOUBLE, myRank^shift, DEFAULT_TAG, MPI_COMM_WORLD);
        MPI_Recv(&rResult, 1, MPI_DOUBLE, myRank^shift, DEFAULT_TAG, MPI_COMM_WORLD, &status);

        myResult += rResult;

        // printf("I am process %d and I have received from process %d\n", myRank, myRank^shift);
    }

    // make some barriers and measure the time
    if (myRank == root){
        printf("Result: %lf\n", myResult);
    }

    // detach communication buffer
    MPI_Buffer_detach(buffer, &bufferSize);

    // finalizing MPI interface
    MPI_Finalize();
}
