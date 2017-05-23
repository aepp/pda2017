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
//
//void nonBlockingCommunication(int* myRandomInts, int myRank, int numProc, int sizeOfRandArray){
//    MPI_Request sendReqO,                       // odd send request handle
//                sendReqE,                       // even send request handle
//                recvReqO,                       // odd recv request handle
//                recvReqE;                       // even recv request handle
//    MPI_Status  statusO,                        // odd communication status
//                statusE;                        // even communication status
//    int myNewRandomInts[sizeOfRandArray],          // elements current process keep
//        randomIntsToSend[sizeOfRandArray],       // elements current process sends to its communication partner
//        phase,                                  // phase (even or odd)
//        partnerE,
//        partnerO;
//
//    // calculate communication partner rank for even and odd phase
//    if(myRank % 2 == 0){
//        partnerE = myRank + 1;
//        partnerO = myRank - 1;
//        if(partnerE >= numProc){
//            partnerE = -1;
//        }
//    } else {
//        partnerE = myRank - 1;
//        partnerO = myRank + 1;
//        if(partnerO >= numProc){
//            partnerO = -1;
//        }
//    }
//
//    // odd-even transportation sort
//    for (phase = 0; phase < numProc; phase++){
//        if (phase % 2 == 0 && partnerE > 0){ // even phase
//            MPI_Isend(myRandomInts, sizeOfRandArray, MPI_INT, partnerE, 1, MPI_COMM_WORLD, &sendReqE);
//
//            MPI_Irecv(myNewRandomInts, sizeOfRandArray, MPI_INT, partnerE, 1, MPI_COMM_WORLD, &recvReqE);
//        } else if (phase % 2 != 0 && partnerO > 0){ // odd phase
//            MPI_Irecv(theirRandomInts, sizeOfRandArray, MPI_INT, partnerO, 1, MPI_COMM_WORLD, &recvReqO);
//
//            MPI_Wait(&recvReqO, &statusO);
//
//            sortArrays(randomInts, theirRandomInts, myRandomInts, sizeOfRandArray);
//
//            MPI_Isend(rand_ints_send, n, MPI_INT, partner, 1, MPI_COMM_WORLD, &sendReqO);
//            MPI_Wait(&sendReqO, &statusO);
//        }
//    }
//}
//
//void persistentCommunication(int* randomInts, int myRank, int numProc, int sizeOfRandArray){
//
//}
