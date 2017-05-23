#include <stdio.h>                          // import of the definitions of the C IO library
#include <string.h>                         // import of the definitions of the string operations
#include <unistd.h>                         // standard unix io library definitions and declarations
#include <errno.h>                          // system error numbers
#include <stdlib.h>                         // for random()
#include <time.h>                           // to seed random generator
#include <math.h>                           // for math functions
#include "mpi.h"                            // import of the MPI definitions

#include "util.h"                           // include own header file

double f1(double x){
    return 1 / (sqrt(2 * x + 1));
}

double f2(double x){
    return log(x);
}

double trapezoidalRuleF1(double a, double b, int n){
    int k;
    double tempSum = 0,
           s = (b - a)/n;

    for(k = 1; k < n; k++) {
        tempSum += f1(a + k * s);
    }
    return (s/2) * (f1(a) + f1(b) + 2 * tempSum);
}

double trapezoidalRuleF2(double a, double b, int n){
    int k;
    double tempSum = 0,
           s = (b - a)/n;

    for(k = 1; k < n; k++) {
        tempSum += f2(a + k * s);
    }
    return ((s/n)/2) * (f2(a) + f2(b) + 2 * tempSum);
}

void fillWithRandomInt(int* array, int size, int maxRandom, int rank)
{
    srandom((unsigned)time(NULL) * rank);
    int i;
    for(i = 0; i < size; i++) {
        // assign random int to each array position
        array[i] = generateRandomInt(maxRandom);
    }
}

int generateRandomInt(int max){
    return random() % max;
}

int cmpFunc (const void *a, const void *b)
{
   return (*(int*)a - *(int*)b);
}

void nonBlockingCommunication(int *myRandomInts, int myRank, int numProc, int sizeOfRandArray){
    MPI_Request sendReqO,                   // odd send request handle
                sendReqE,                   // even send request handle
                recvReqO,                   // odd recv request handle
                recvReqE;                   // even recv request handle
    MPI_Status statusO,                     // odd communication status
               statusE;                     // even communication status
    int theirRandomInts[sizeOfRandArray],   // partners random array
        phase,                              // phase (even or odd)
        partnerE,                           // communication partner for the even phase
        partnerO;                           // communication partner for the odd phase

    // calculate communication partner rank for even and odd phase
    if(myRank % 2 == 0){
        partnerE = myRank + 1;
        partnerO = myRank - 1;
        if(partnerE >= numProc){
            partnerE = -1;
        }
    } else {
        partnerE = myRank - 1;
        partnerO = myRank + 1;
        if(partnerO >= numProc){
            partnerO = -1;
        }
    }

    // odd-even transportation sort
    for (phase = 0; phase < numProc; phase++){
        if (phase % 2 == 0 && partnerE >= 0){ // even phase
            MPI_Isend(myRandomInts, sizeOfRandArray, MPI_INT, partnerE, 1, MPI_COMM_WORLD, &sendReqE);
            MPI_Irecv(theirRandomInts, sizeOfRandArray, MPI_INT, partnerE, 1, MPI_COMM_WORLD, &recvReqE);
            // wait until receive finished
            MPI_Wait(&recvReqE, &statusE);
            // merge own and received arrays, sort and decide which half to keep
            getMyNewArray(myRandomInts, theirRandomInts, sizeOfRandArray, partnerE, myRank);
        } else if (phase % 2 != 0 && partnerO >= 0){ // odd phase
            MPI_Isend(myRandomInts, sizeOfRandArray, MPI_INT, partnerO, 1, MPI_COMM_WORLD, &sendReqO);
            MPI_Irecv(theirRandomInts, sizeOfRandArray, MPI_INT, partnerO, 1, MPI_COMM_WORLD, &recvReqO);
            // wait until send finished
            MPI_Wait(&recvReqO, &statusO);
            // merge own and received arrays, sort and decide which half to keep
            getMyNewArray(myRandomInts, theirRandomInts, sizeOfRandArray, partnerO, myRank);
        }
    }
}

void persistentCommunication(int* myRandomInts, int myRank, int numProc, int sizeOfRandArray){
    MPI_Request sendReqO,                   // odd send request handle
                sendReqE,                   // even send request handle
                recvReqO,                   // odd recv request handle
                recvReqE;                   // even recv request handle
    MPI_Status statusO,                     // odd communication status
               statusE;                     // even communication status
    int theirRandomInts[sizeOfRandArray],   // partners random array
        phase,                              // phase (even or odd)
        partnerE,                           // communication partner for the even phase
        partnerO;                           // communication partner for the odd phase

    // calculate communication partner rank for even and odd phase
    if(myRank % 2 == 0){
        partnerE = myRank + 1;
        partnerO = myRank - 1;
        if(partnerE >= numProc){
            partnerE = -1;
        }
    } else {
        partnerE = myRank - 1;
        partnerO = myRank + 1;
        if(partnerO >= numProc){
            partnerO = -1;
        }
    }

    // initialize persistent communication
    MPI_Send_init(myRandomInts, sizeOfRandArray, MPI_INT, partnerE, 1, MPI_COMM_WORLD, &sendReqE);
    MPI_Recv_init(theirRandomInts, sizeOfRandArray, MPI_INT, partnerE, 1, MPI_COMM_WORLD, &recvReqE);
    MPI_Send_init(myRandomInts, sizeOfRandArray, MPI_INT, partnerO, 1, MPI_COMM_WORLD, &sendReqO);
    MPI_Recv_init(theirRandomInts, sizeOfRandArray, MPI_INT, partnerO, 1, MPI_COMM_WORLD, &recvReqO);

    // odd-even transportation sort
    for (phase = 0; phase < numProc; phase++){
        if (phase % 2 == 0 && partnerE >= 0){ // even phase
            // start communication
            MPI_Start(&sendReqE);
            MPI_Start(&recvReqE);
            // wait until requests are finished
            MPI_Wait(&sendReqE, &statusE);
            MPI_Wait(&recvReqE, &statusE);

            // merge own and received arrays, sort and decide which half to keep
            getMyNewArray(myRandomInts, theirRandomInts, sizeOfRandArray, partnerE, myRank);
        } else if (phase % 2 != 0 && partnerO >= 0){ // odd phase
            // start communication
            MPI_Start(&sendReqO);
            MPI_Start(&recvReqO);
            // wait until requests are finished
            MPI_Wait(&sendReqO, &statusO);
            MPI_Wait(&recvReqO, &statusO);
            // merge own and received arrays, sort and decide which half to keep
            getMyNewArray(myRandomInts, theirRandomInts, sizeOfRandArray, partnerO, myRank);
        }
    }

    MPI_Request_free(&sendReqO);
    MPI_Request_free(&sendReqE);
    MPI_Request_free(&recvReqO);
    MPI_Request_free(&recvReqE);
}

void getMyNewArray(int *myRandomInts, int *theirRandomInts, int sizeOfRandArray, int partner, int myRank){
    int mergedArray[2 * sizeOfRandArray];
    // merge own random array
    memcpy(mergedArray, myRandomInts, sizeOfRandArray * sizeof(int));
    // and partners random array into one
    memcpy(mergedArray + sizeOfRandArray, theirRandomInts, sizeOfRandArray * sizeof(int));
    // sort the merged array
    qsort(mergedArray, 2 * sizeOfRandArray, sizeof(int), cmpFunc);

    // apply the exchange rule
    if (myRank < partner){ // if I'm left from my partner
        // keep the left part of the merged sorted array (smaller values)
        memcpy(myRandomInts, mergedArray, sizeOfRandArray * sizeof(int));
    } else { // if I'm the right neighbour of my partner
        // keep the right part of the merged sorted array (bigger values)
        memcpy(myRandomInts, mergedArray + sizeOfRandArray, sizeOfRandArray * sizeof(int));
    }
}