/************************************************************************/
/* Author: Aleksandr Epp <aleksandr.epp@gmail.com>                      */
/* matriclenumber: 6002853                                              */
/* Assignment : 1                                                       */
/* Task : 2    a)                                                       */
/*                                                                      */
/* Description:                                                         */
/*                                                                      */
/*                                                                      */
/************************************************************************/

#include "mpi.h"         // import of the MPI definitions
#include <stdio.h>         // import of the definitions of the C IO library
#include <string.h>     // import of the definitions of the string operations
#include <unistd.h>        // standard unix io library definitions and declarations
#include <errno.h>        // system error numbers

#include "task2a.h"        // include own header file
#include "util.h"         // include assignment utils

//#define MPI_WTIME_IS_GLOBAL 1

void task2a(int argc, char* argv[ ], int maxRandom, int arraySize)
{
    // rank of the process
    int myRank;
    // rank of the process sent the message
    int theirRank;

    // messages tags
    int TAG_MIN = 1;    // tag for messages with min
    int TAG_MAX = 2;    // tag for messages with max
    int TAG_SUM = 3;    // tag for messages with sum
    
    // number of proccesses envolved
    int numProc; 
    // iterator variable
    int i;

    double timeStart, timeStop; // begin and end times to measure execution time

    // initializing of MPI-Interface
    MPI_Init(&argc, &argv);

    MPI_Request request; //capture request of a MPI_Send
    MPI_Status status; //capture status of a MPI_Send
    
    //get your rank
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // get the number of processes running
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);

    // start time measurement when all initialization is done and calculation begins
    timeStart = MPI_Wtime();

    // array to operate on
    int randomInts[arraySize];
    // min, max and sum of the array calculated by current process
    int myMin, myMax, mySum;
    // min, max and sum of the array recieved form another process
    int rMin = 42, rMax, rSum;
    
    // fill with random integers
    fillWithRandomInt(randomInts, arraySize, maxRandom);
    
    // print initialized array
//    for(i = 0; i < arraySize; i++) {
//        printf("%d ", randomInts[i]);
//        printf("\n");
//    }
    
    // calculate min, max and sum of the own array
    myMin = getMin(randomInts, arraySize);
    myMax = getMax(randomInts, arraySize);
    mySum = getSum(randomInts, arraySize);

    
    if(myRank == 0){ // if I'm the first process:
        // send my results to next process
        MPI_Send(&myMin, 1, MPI_INT, myRank + 1, TAG_MIN, MPI_COMM_WORLD);
        MPI_Send(&myMax, 1, MPI_INT, myRank + 1, TAG_MAX, MPI_COMM_WORLD);
        MPI_Send(&mySum, 1, MPI_INT, myRank + 1, TAG_SUM, MPI_COMM_WORLD);
        
//        printf("I'm %d. My Min=%d. I sent Min=%d to %d\n", myRank, myMin, myMin, myRank + 1);
        
        // wait for results from the last process
        MPI_Recv(&rMin, 1, MPI_INT, numProc - 1, TAG_MIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&rMax, 1, MPI_INT, numProc - 1, TAG_MAX, MPI_COMM_WORLD, &status);
        MPI_Recv(&rSum, 1, MPI_INT, numProc - 1, TAG_SUM, MPI_COMM_WORLD, &status);

        // end time measurement when first process receives the results from the last one
        timeStop = MPI_Wtime();
        
//        printf("I'm %d. I recieved Min=%d from %d\n", myRank, rMin, numProc - 1);
        // and print the overall result
        printf("Overall results:\n\n");
        printf("Min: %d\n", rMin);
        printf("Max: %d\n", rMax);
        printf("Sum: %d\n", rSum);
        printf("\nOverall execution time: %f\n", timeStop - timeStart);
        
    } else if(myRank < numProc - 1) { // if I'm process not first and not las process:
        // recieve results from the previous process
        MPI_Recv(&rMin, 1, MPI_INT, myRank - 1, TAG_MIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&rMax, 1, MPI_INT, myRank - 1, TAG_MAX, MPI_COMM_WORLD, &status);
        MPI_Recv(&rSum, 1, MPI_INT, myRank - 1, TAG_SUM, MPI_COMM_WORLD, &status);
        
        //MPI_Wait(&request, &status);
        
//        printf("I'm %d. My Min=%d. I recieved Min=%d from %d\n", myRank, myMin, rMin, myRank - 1);
        
        // calculate min, max and sum considering recieved results
        if(rMin < myMin) {
            myMin = rMin;
        }
        if(rMax > myMax) {
            myMax = rMax;
        }
        mySum += rSum;
        
        // and send my results to next process
        MPI_Send(&myMin, 1, MPI_INT, myRank + 1, TAG_MIN, MPI_COMM_WORLD);
        MPI_Send(&myMax, 1, MPI_INT, myRank + 1, TAG_MAX, MPI_COMM_WORLD);
        MPI_Send(&mySum, 1, MPI_INT, myRank + 1, TAG_SUM, MPI_COMM_WORLD);
        
//        printf("I'm %d. I sent Min=%d to %d\n", myRank, myMin, myRank + 1);
    } else { // if I'm the last process:
        // recieve results from the previous process
        MPI_Recv(&rMin, 1, MPI_INT, myRank - 1, TAG_MIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&rMax, 1, MPI_INT, myRank - 1, TAG_MAX, MPI_COMM_WORLD, &status);
        MPI_Recv(&rSum, 1, MPI_INT, myRank - 1, TAG_SUM, MPI_COMM_WORLD, &status);

//        printf("I'm %d. My Min=%d. I recieved Min=%d from %d\n", myRank, myMin, rMin, myRank - 1);
        
        // calculate min, max and sum considering recieved results
        if(rMin < myMin) {
            myMin = rMin;
        }
        if(rMax > myMax) {
            myMax = rMax;
        }
        mySum += rSum;
        
        // and send my results to the first process
        MPI_Send(&myMin, 1, MPI_INT, 0, TAG_MIN, MPI_COMM_WORLD);
        MPI_Send(&myMax, 1, MPI_INT, 0, TAG_MAX, MPI_COMM_WORLD);
        MPI_Send(&mySum, 1, MPI_INT, 0, TAG_SUM, MPI_COMM_WORLD);
        
//        printf("I'm %d. I sent Min=%d to %d\n", myRank, myMin, 0);
    }

    // finalizing MPI interface
    MPI_Finalize();    
}
