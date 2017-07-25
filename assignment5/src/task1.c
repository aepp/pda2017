/************************************************************************/
/* Author: Aleksandr Epp <aleksandr.epp@gmail.com>                      */
/* Matriclenumber: 6002853                                              */
/* Assignment : 5                                                       */
/* Task : 1                                                             */
/*                                                                      */
/* Description:                                                         */
/*                                                                      */
/* Distributed image processing filter                                  */
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

#define TAG1    1

void task1(int argc, char* argv[], char* inputImgFilePath, int filterType, int filterStrength)
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
    MPI_Offset totalPixelCount;                 // total pixel count in the input image
    MPI_Status readWriteStatus;

    // initializing of MPI-Interface
    MPI_Init(&argc, &argv);

    MPI_Status status; //capture status of a MPI_Send

    // get your rank
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // get the number of processes running
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);

    error = MPI_File_open(MPI_COMM_WORLD, inputImgFilePath, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

    // get amount of pixels in the image
    MPI_File_get_size(fh, &totalPixelCount);

    // determine how many rows each process reads
    rowsPerProcessCount = totalPixelCount/(1280 * numProc * sizeof(unsigned char));

    // determined amount of rows, 1280 pixels in each row
    unsigned char myImgRows[rowsPerProcessCount][1280];
    // image rows with 2 pixel padding on the top and bottom sides
    unsigned char myTBPaddedImgRows[rowsPerProcessCount + 4][1280];
    // image rows with 2 pixel padding on each side
    unsigned char myCompletelyPaddedImgRows[rowsPerProcessCount + 4][1280 + 4];
    // store filtered result here
    unsigned char myFilteredImgRows[rowsPerProcessCount][1280];

    // last row from previous process to prepend to own rows
    unsigned char lastRowFormPrevProc[2][1280];
    // first row from next process to append to own rows
    unsigned char firstRowFormNextProc[2][1280];

    MPI_File_read_ordered(fh, myImgRows, rowsPerProcessCount * 1280, MPI_UNSIGNED_CHAR, &readWriteStatus);
    MPI_File_close(&fh);

    // iterator for filter strength
    int m;
    for(m = 0; m < filterStrength; m++){
        if(numProc > 1){
            if(myRank != 0 && myRank != numProc - 1){
                // send my last row to the next process
                MPI_Send(myImgRows[rowsPerProcessCount - 2], 2 * 1280, MPI_UNSIGNED_CHAR, myRank + 1, TAG1, MPI_COMM_WORLD);
                // send my first row to the previous process
                MPI_Send(myImgRows[0], 2 * 1280, MPI_UNSIGNED_CHAR, myRank - 1, TAG1, MPI_COMM_WORLD);
                // receive one row from process before
                MPI_Recv(lastRowFormPrevProc, 2 * 1280, MPI_UNSIGNED_CHAR, myRank - 1, TAG1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // receive one row from process after
                MPI_Recv(firstRowFormNextProc, 2 * 1280, MPI_UNSIGNED_CHAR, myRank + 1, TAG1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                // prepend received rows to my padded rows
                memcpy(myTBPaddedImgRows, lastRowFormPrevProc, 2 * 1280 * sizeof(unsigned char));

                // append received rows to my padded rows
                memcpy((unsigned char*)((uintptr_t)myTBPaddedImgRows + (2 + rowsPerProcessCount) * 1280), firstRowFormNextProc, 2 * 1280 * sizeof(unsigned char));
            } else if(myRank == 0){
                // send my last row to the next process
                MPI_Send(myImgRows[rowsPerProcessCount - 2], 2 * 1280, MPI_UNSIGNED_CHAR, myRank + 1, TAG1, MPI_COMM_WORLD);
                // receive one row from process after
                MPI_Recv(firstRowFormNextProc, 2 * 1280, MPI_UNSIGNED_CHAR, myRank + 1, TAG1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // prepend empty (black) rows to my padded rows
                memcpy(myTBPaddedImgRows, calloc(2 * 1280, sizeof(unsigned char)), 2 * 1280 * sizeof(unsigned char));
                // append received rows to my padded rows
                memcpy((unsigned char*)((uintptr_t)myTBPaddedImgRows + (2 + rowsPerProcessCount) * 1280), firstRowFormNextProc, 2 * 1280 * sizeof(unsigned char));
            } else if(myRank == numProc - 1){
                // send my first row to the previous process
                MPI_Send(myImgRows[0], 2 * 1280, MPI_UNSIGNED_CHAR, myRank - 1, TAG1, MPI_COMM_WORLD);
                // receive one row from process before
                MPI_Recv(lastRowFormPrevProc, 2 * 1280, MPI_UNSIGNED_CHAR, myRank - 1, TAG1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // prepend received rows to my padded rows
                memcpy(myTBPaddedImgRows, lastRowFormPrevProc, 2 * 1280 * sizeof(unsigned char));
                // append  empty (black) rows to my padded rows
                memcpy((unsigned char*)((uintptr_t)myTBPaddedImgRows + (2 + rowsPerProcessCount) * 1280), calloc(2 * 1280, sizeof(unsigned char)), 2 * 1280 * sizeof(unsigned char));
            }
        } else {
            // prepend empty (black) rows to my padded rows
            memcpy(myTBPaddedImgRows, calloc(2 * 1280, sizeof(unsigned char)), 2 * 1280 * sizeof(unsigned char));
            // append  empty (black) rows to my padded rows
            memcpy((unsigned char*)((uintptr_t)myTBPaddedImgRows + (2 + rowsPerProcessCount) * 1280), calloc(2 * 1280, sizeof(unsigned char)), 2 * 1280 * sizeof(unsigned char));
        }
        // copy own image part (rows) to padded rows array,
        // where 2 rows from previous process are prepended
        // and 2 rows from next process are appended
        memcpy((unsigned char*)((uintptr_t)myTBPaddedImgRows + 2 * 1280), myImgRows, rowsPerProcessCount * 1280 * sizeof(unsigned char));

        // check if top and bottom padding applied properly
//        error = MPI_File_open(MPI_COMM_WORLD, "./TBPaddedImg.gray", MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
//        MPI_File_write_ordered(fh, myTBPaddedImgRows, (rowsPerProcessCount + 4) * 1280, MPI_UNSIGNED_CHAR, &readWriteStatus);
//        MPI_File_close(&fh);

        addLeftRightPadding(myTBPaddedImgRows, myCompletelyPaddedImgRows, rowsPerProcessCount + 4);

        // check left and right padding applied properly
//        error = MPI_File_open(MPI_COMM_WORLD, "./CompletelyPaddedImgRows.gray", MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
//        MPI_File_write_ordered(fh, myCompletelyPaddedImgRows, (rowsPerProcessCount + 4) * (1280 + 4), MPI_UNSIGNED_CHAR, &readWriteStatus);
//        MPI_File_close(&fh);

        applyFilter(myCompletelyPaddedImgRows, myFilteredImgRows, filterType, rowsPerProcessCount);
        memcpy(myImgRows, myFilteredImgRows, rowsPerProcessCount * 1280 * sizeof(unsigned char));
    }

    // save filtered image
    error = MPI_File_open(MPI_COMM_WORLD, "./filtered.gray", MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_File_write_ordered(fh, myFilteredImgRows, rowsPerProcessCount * 1280, MPI_UNSIGNED_CHAR, &readWriteStatus);
    MPI_File_close(&fh);
    // finalizing MPI interface
    MPI_Finalize();
}
