/************************************************************************/
/* Author: Aleksandr Epp <aleksandr.epp@gmail.com>                      */
/* Matriclenumber: 6002853                                              */
/* Assignment : 4                                                       */
/* Task : 1                                                             */
/*                                                                      */
/* Description:                                                         */
/*                                                                      */
/* Shared Sort on 2D Mesh                                               */
/*                                                                      */
/************************************************************************/

#include "mpi.h"                                // import of the MPI definitions
#include <stdio.h>                              // import of the definitions of the C IO library
#include <string.h>                             // import of the definitions of the string operations
#include <stdlib.h>                             // to dynamically allocate array
#include <unistd.h>                             // standard unix io library definitions and declarations
#include <errno.h>                              // system error numbers
#include <math.h>                               // for log()

#include "task1.h"                              // include own header file
#include "util.h"                               // include assignment utils

#define MAX_RANDOM      100                     // max random
#define MESH_DIMS        2                       // amount of mesh dimensions (2D Mesh)

// topology constants
#define NOWRAP          0                       // non-cyclic
#define REORDER         1                       // reordering (allowed = 1, not allowed = 0)
void task1(int argc, char* argv[], int randomIntCount)
{
    int myRank,                                 // rank of the process
        theirRank,                              // rank of the process sent the message
        DEFAULT_TAG = 1,                        // tag for messages with min
        numProc,                                // number of processes
        i,                                      // iterator variable
        root = 0,                               // root process
        meshDimSize;                            // single dimension of a squared grid


    // initializing of MPI-Interface
    MPI_Init(&argc, &argv);

    //capture status of a MPI_Send
    MPI_Status status;

    // get your rank
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // get the number of processes running
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);

    MPI_Comm gridComm, // communicator associated with grid
             gridColComm, // communicator associated with grid columns
             gridRowComm; // communicator associated with grid rows

    // calculate the size of dimensions depending on processes involved
    meshDimSize = (int)sqrt(numProc);

    int dimSize[MESH_DIMS] = {meshDimSize, meshDimSize}, // define sizes for each dimension
        wrap[MESH_DIMS] = {NOWRAP, NOWRAP}, // non-cycling dimensions
        remainDimsCols[MESH_DIMS] = {1, 0}, // 1st dim varies, 2nd dim fixed
        remainDimsRows[MESH_DIMS] = {0, 1}; // 1st dim fixed, 2nd dim varies

    // create communication grid
    MPI_Cart_create(
        MPI_COMM_WORLD,
        MESH_DIMS,
        dimSize,
        wrap,
        REORDER,
        &gridComm
    );

    // split grid in columns
    MPI_Cart_sub(gridComm, remainDimsCols, &gridColComm);

    // split grid in rows
    MPI_Cart_sub(gridComm, remainDimsRows, &gridRowComm);

    int phase;                      // phase counter
    int gridColRank,                // process rank in grid column communicator
        gridRowRank,                // process rank in grid row communicator
        gridColCommSize,            // amount of processes in grid column communicator
        gridRowCommSize,            // amount of processes in grid row communicator
        randomInts[randomIntCount], // random array
        sortOrder = 0;              // default sort order is descending for odd rows

    // get your local rank
    MPI_Comm_rank(gridColComm, &gridColRank);
    MPI_Comm_rank(gridRowComm, &gridRowRank);

    // get the number of processes in comm
    MPI_Comm_size(gridColComm, &gridColCommSize);
    MPI_Comm_size(gridRowComm, &gridRowCommSize);

    // fill with random integers
    fillWithRandomInt(randomInts, randomIntCount, MAX_RANDOM, myRank);

    // decide how to sort depending on which row the process belong to
    if(gridColRank % 2 == 0){ // even
//    if((myRank / meshDimSize) % 2 == 0) {
        sortOrder = 1;
    }

    // shared sort
    for(phase = 0; phase <= 2 * log2(meshDimSize) + 1; phase++){
        // sort rows
        oddEvenTranspositionSort(randomInts, gridRowComm, randomIntCount, sortOrder);
        // sort columns
        oddEvenTranspositionSort(randomInts, gridColComm, randomIntCount, 1);
    }

    // use task 1 from assignment 1 to make ordered output
    int theirRandomInts[randomIntCount], j;
    if(myRank == root){ // if I'm process with rank 0
        printf("me: %d -> [", myRank);
        for(j = 0; j < randomIntCount - 1; j++){
            printf("%d, ", randomInts[j]);
        }
        printf("%d]\n", randomInts[randomIntCount - 1]);
        // receive messages from all other processes
        for(i = 1; i < numProc; i++){
            // receive the process rank
            MPI_Recv(theirRandomInts, randomIntCount, MPI_INT, i, DEFAULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // print the received message to console
            printf("me: %d -> [", i);
            for(j = 0; j < randomIntCount - 1; j++){
                printf("%d, ", theirRandomInts[j]);
            }
            printf("%d]\n", theirRandomInts[randomIntCount - 1]);
        }
    } else { // if I'm process with rank > 0, send message to process with rank 0
        // send my rank
        MPI_Send(randomInts, randomIntCount, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD);
    }

    // finalizing MPI interface
    MPI_Finalize();
}
