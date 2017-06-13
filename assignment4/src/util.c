#include <stdio.h>                          // import of the definitions of the C IO library
#include <string.h>                         // import of the definitions of the string operations
#include <unistd.h>                         // standard unix io library definitions and declarations
#include <errno.h>                          // system error numbers
#include <stdlib.h>                         // for random()
#include <time.h>                           // to seed random generator
#include <math.h>                           // for math functions
#include "mpi.h"                            // import of the MPI definitions

#include "util.h"                           // include own header file

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

int cmpFuncASC (const void *a, const void *b)
{
   return (*(int*)a - *(int*)b);
}

int cmpFuncDESC (const void *a, const void *b)
{
   return (*(int*)b - *(int*)a);
}


void oddEvenTranspositionSort(int *myRandomInts, MPI_Comm comm, int sizeOfRandArray, int sortOrder){
    int myRank,                             // process rank in communicator
        numProc,                            // amount of processes in communicator
        theirRandomInts[sizeOfRandArray],   // partners random array
        phase,                              // phase (even or odd)
        partnerE,                           // communication partner for the even phase
        partnerO;                           // communication partner for the odd phase

    // get your local rank
    MPI_Comm_rank(comm, &myRank);

    // get the number of processes in comm
    MPI_Comm_size(comm, &numProc);

    MPI_Request sendReqO,                   // odd send request handle
                sendReqE,                   // even send request handle
                recvReqO,                   // odd recv request handle
                recvReqE;                   // even recv request handle
    MPI_Status statusO,                     // odd communication status
               statusE;                     // even communication status

    // decide how to sort
    if(sortOrder == 1){
        qsort(myRandomInts, sizeOfRandArray, sizeof(int), cmpFuncASC);
    } else {
        qsort(myRandomInts, sizeOfRandArray, sizeof(int), cmpFuncDESC);
    }

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
            MPI_Isend(myRandomInts, sizeOfRandArray, MPI_INT, partnerE, 1, comm, &sendReqE);
            MPI_Irecv(theirRandomInts, sizeOfRandArray, MPI_INT, partnerE, 1, comm, &recvReqE);
            // wait until receive finished
            MPI_Wait(&recvReqE, &statusE);
            // merge own and received arrays, sort and decide which half to keep
            getMyNewArray(myRandomInts, theirRandomInts, sizeOfRandArray, partnerE, myRank, sortOrder);
        } else if (phase % 2 != 0 && partnerO >= 0){ // odd phase
            MPI_Isend(myRandomInts, sizeOfRandArray, MPI_INT, partnerO, 1, comm, &sendReqO);
            MPI_Irecv(theirRandomInts, sizeOfRandArray, MPI_INT, partnerO, 1, comm, &recvReqO);
            // wait until send finished
            MPI_Wait(&recvReqO, &statusO);
            // merge own and received arrays, sort and decide which half to keep
            getMyNewArray(myRandomInts, theirRandomInts, sizeOfRandArray, partnerO, myRank, sortOrder);
        }
    }
}

void getMyNewArray(int *myRandomInts, int *theirRandomInts, int sizeOfRandArray, int partner, int myRank, int sortOrder){
    int mergedArray[2 * sizeOfRandArray];
    // merge own random array
    memcpy(mergedArray, myRandomInts, sizeOfRandArray * sizeof(int));
    // and partners random array into one
    memcpy(mergedArray + sizeOfRandArray, theirRandomInts, sizeOfRandArray * sizeof(int));

    // sort the merged array
    // decide how to sort
    if(sortOrder == 1){
        qsort(mergedArray, 2 * sizeOfRandArray, sizeof(int), cmpFuncASC);
    } else {
        qsort(mergedArray, 2 * sizeOfRandArray, sizeof(int), cmpFuncDESC);
    }

    // apply the exchange rule
    if (myRank < partner){ // if I'm left from my partner
        // keep the left part of the merged sorted array (smaller values)
        memcpy(myRandomInts, mergedArray, sizeOfRandArray * sizeof(int));
    } else { // if I'm the right neighbour of my partner
        // keep the right part of the merged sorted array (bigger values)
        memcpy(myRandomInts, mergedArray + sizeOfRandArray, sizeOfRandArray * sizeof(int));
    }
}