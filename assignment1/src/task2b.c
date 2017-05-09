/************************************************************************/
/* Author: Aleksandr Epp <aleksandr.epp@gmail.com>                      */
/* Matriclenumber: 6002853                                              */
/* Assignment : 1                                                       */
/* Task : 2 b)                                                          */
/*                                                                      */
/* Description:                                                         */
/*                                                                      */
/* Same calculations as in Task 2 a) but processes communicate in a     */
/* tree way.                                                            */
/*                                                                      */
/************************************************************************/

#include "mpi.h"        // import of the MPI definitions
#include <stdio.h>      // import of the definitions of the C IO library
#include <string.h>     // import of the definitions of the string operations
#include <unistd.h>     // standard unix io library definitions and declarations
#include <errno.h>      // system error numbers
#include <math.h>       // for log2()
#include "task2b.h"     // include own header file
#include "util.h"       // include assignment utils

void task2b(int argc, char* argv[ ], int maxRandom, int arraySize)
{
    // rank of the process
    int myRank;
    // rank of the process sent the message
    int theirRank;

    // messages tags
    int TAG_MIN = 1;    // tag for messages with min
    int TAG_MAX = 2;    // tag for messages with max
    int TAG_SUM = 3;    // tag for messages with sum

    // number of processes
    int numProc;
    // iterator variable
    int i;

    double timeStart, timeStop; // begin and end times to measure execution time

    // array to operate on
    int randomInts[arraySize],
        // min, max and sum of the array calculated by current process
        myMin, myMax, mySum,
        // min, max and sum of the array received form another process
        rMin, rMax, rSum,
        // amount of communication levels
        levelCount,
        // level of current process, set to 0 for the first process
        myLevel = 0;

    // initializing of MPI-Interface
    MPI_Init(&argc, &argv);

    // start time measurement when initialization
    timeStart = MPI_Wtime();

    // start time measurement when all initialization is done and calculation begins
    timeStart = MPI_Wtime();

    MPI_Status status; //capture status of a MPI_Send

    //get your rank
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // get the number of processes running
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);

    // determine the amount of communication levels
    levelCount = (int)log2l(numProc);

    // fill with random integers
    fillWithRandomInt(randomInts, arraySize, maxRandom, myRank);

    // calculate min, max and sum of the own array
    myMin = getMin(randomInts, arraySize);
    myMax = getMax(randomInts, arraySize);
    mySum = getSum(randomInts, arraySize);

//    for (i = 0; i < arraySize; i++) {
//        printf("Me %d: array[%d]=%d\n", myRank, i, randomInts[i]);
//    }

//    printf("Me %d: myMax=%d\n", myRank, myMax);

    for(i = 0; i < levelCount; i++) {
        // determine level for each process
        // level of process 0 will stay 0
        if(myRank >= (int)powl(2, i) && myRank < pow(2, i + 1)){
            myLevel = i + 1;
        }
    }

//    printf("Me %d! My level is %d\n", myRank, myLevel);

    if(myRank >= numProc / 2) { // all processes with rank > 0 send
//        printf("Me %d! Send to %d\n", myRank, myRank - (int)powl(2, myLevel - 1));
        MPI_Send(&myMin, 1, MPI_INT, myRank - (int)powl(2, myLevel - 1), TAG_MIN, MPI_COMM_WORLD);
        MPI_Send(&myMax, 1, MPI_INT, myRank - (int)powl(2, myLevel - 1), TAG_MAX, MPI_COMM_WORLD);
        MPI_Send(&mySum, 1, MPI_INT, myRank - (int)powl(2, myLevel - 1), TAG_SUM, MPI_COMM_WORLD);
    } else if(myRank < numProc / 2 && myRank != 0) { // processes with rank < numProc / 2 receive and send
        for(i = myLevel; i < levelCount; i++) {
//            printf("Me %d! Receive from %d\n", myRank, myRank + (int)powl(2, i));
            MPI_Recv(&rMin, 1, MPI_INT, myRank + (int)powl(2, i), TAG_MIN, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&rMax, 1, MPI_INT, myRank + (int)powl(2, i), TAG_MAX, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&rSum, 1, MPI_INT, myRank + (int)powl(2, i), TAG_SUM, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            if(rMin < myMin) {
                myMin = rMin;
            }
            if(rMax > myMax) {
                myMax = rMax;
            }
            mySum += rSum;
        }

//        printf("Me %d! Send to %d\n", myRank, myRank - (int)powl(2, myLevel - 1));
        MPI_Send(&myMin, 1, MPI_INT, myRank - (int)powl(2, myLevel - 1), TAG_MIN, MPI_COMM_WORLD);
        MPI_Send(&myMax, 1, MPI_INT, myRank - (int)powl(2, myLevel - 1), TAG_MAX, MPI_COMM_WORLD);
        MPI_Send(&mySum, 1, MPI_INT, myRank - (int)powl(2, myLevel - 1), TAG_SUM, MPI_COMM_WORLD);
    } else {
         for(i = myLevel; i < levelCount; i++) {
//            printf("Me %d! Receive from %d\n", myRank, myRank + (int)powl(2, i));
            MPI_Recv(&rMin, 1, MPI_INT, myRank + (int)powl(2, i), TAG_MIN, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&rMax, 1, MPI_INT, myRank + (int)powl(2, i), TAG_MAX, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&rSum, 1, MPI_INT, myRank + (int)powl(2, i), TAG_SUM, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//            printf("Me %d! Received max=%d from %d\n", myRank, rMax, myRank + (int)powl(2, i));

            if(rMin < myMin) {
                myMin = rMin;
            }
            if(rMax > myMax) {
                myMax = rMax;
            }
            mySum += rSum;
        }

        // end time measurement when first process receives the results
        timeStop = MPI_Wtime();

        // and print the overall result
        printf("Overall results:\n\n");
        printf("Min: %d\n", myMin);
        printf("Max: %d\n", myMax);
        printf("Sum: %d\n", mySum);
        printf("\nOverall execution time: %f\n", timeStop - timeStart);
    }

    // finalizing MPI interface
    MPI_Finalize();
}