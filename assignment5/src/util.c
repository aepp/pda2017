#include <stdio.h>                          // import of the definitions of the C IO library
#include <string.h>                         // import of the definitions of the string operations
#include <unistd.h>                         // standard unix io library definitions and declarations
#include <errno.h>                          // system error numbers
#include <stdlib.h>                         // for random()
#include <time.h>                           // to seed random generator
#include <math.h>                           // for math functions
#include "mpi.h"                            // import of the MPI definitions

#include "util.h"                           // include own header file

void applyFilter(unsigned char input[][1280 + 4], unsigned char output[][1280], int filterType, int rowsPerProcessCount)
{
    double blurFilter[5][5] = {
        {0./37, 0./37, 1./37, 0./37, 0./37},
        {0./37, 2./37, 4./37, 2./37, 0./37},
        {1./37, 4./37, 9./37, 4./37, 1./37},
        {0./37, 2./37, 4./37, 2./37, 0./37},
        {0./37, 0./37, 1./37, 0./37, 0./37}
    };
    double sharpenFilter[5][5] = {
        {0, 0, 0, 0, 0},
        {0, 0, -1, 0, 0},
        {0, -1, 5, -1, 0},
        {0, 0, -1, 0, 0},
        {0, 0, 0, 0, 0}
    };
    double reliefFilter[5][5] = {
        {0, 0, 0, 0, 0},
        {0, -2, -1, 0, 0},
        {0, -1, 1, 1, 0},
        {0, 0, 1, 2, 0},
        {0, 0, 0, 0, 0}
    };
    double edgeFilter[5][5] = {
        {0./4, 0./4, 0./4, 0./4, 0./4},
        {0./4, 1./4, 2./4, 1./4, 0./4},
        {0./4, 2./4, -12./4, 2./4, 0./4},
        {0./4, 1./4, 2./4, 1./4, 0./4},
        {0./4, 0./4, 0./4, 0./4, 0./4},
    };

    switch(filterType){
        case 1: // blur
            getFilterResult(input, output, blurFilter, rowsPerProcessCount);
            break;
        case 2: // sharpen
            getFilterResult(input, output, sharpenFilter, rowsPerProcessCount);
            break;
        case 3: // relief
            getFilterResult(input, output, reliefFilter, rowsPerProcessCount);
            break;
        case 4: // edge
            getFilterResult(input, output, edgeFilter, rowsPerProcessCount);
            break;
        default:
            exit(1);
    }
}

void getFilterResult(unsigned char input[][1280 + 4], unsigned char output[][1280], double filter[][5], int rowsPerProcessCount)
{
    int x, y;   // pixel iterators

    for(y = 2; y < rowsPerProcessCount + 2; y++){
        for(x = 2; x < 1280 + 2; x++){
            output[y - 2][x - 2] = getSinglePixelFilterResult(input, filter, y, x);
        }
    }
}

unsigned char getSinglePixelFilterResult(unsigned char input[][1280 + 4], double filter[][5], int y, int x)
{
    int v, u,   // filter iterators
        k = 2;  // just for fun

    double result = 0.0;

    for(v = 0; v <= 2 * k; v++){
        for(u = 0; u <= 2 * k; u++){
            result += filter[v][u] * (double)(input[y + v - k][x + u - k]);
        }
    }

    if(result < 0.0){
        result = 0;
    } else if(result > 255.0){
        result = 255;
    }
    return (unsigned char)result;
}

void addLeftRightPadding(unsigned char input[][1280], unsigned char paddedInput[][1280 + 4], int rowsPerProcessCountPadded)
{
    int x, y;   // pixel iterators

    for(y = 0; y < rowsPerProcessCountPadded; y++){
        // add left padding
        for(x = 0; x < 2; x++){
            paddedInput[y][x] = 0;
        }

        // copy original image rows to the padded array
        for(x = 2; x < 1280 + 2; x++){
            paddedInput[y][x] = input[y][x - 2];
        }

        // add right padding
        for(x = 1280 + 2; x < 1280 + 4; x++){
            paddedInput[y][x] = 0;
        }
    }
}

void startMultiplications(
    MPI_Comm intercomms[],
    MPI_Request openRequests[],
    MPI_Datatype blockTypes[],
    int workerCounts[],
    int activeMultCount,
    int matrixDimensions[]
)
{
    int i,j,                        // iterators
        totalElementCount,          // total amount of elements in the matrix
        error = 0,
        finishedMultiplications;

    char userInputMatrixAFileName[100],      // file path for Matrix A
         userInputMatrixBFileName[100];      // file path for Matrix B

    MPI_File fh;                    // file handle to open/save files
    MPI_Status status;              // communication status
    MPI_Offset matrixFileSize;      // size of the matrix file

    // read in matrices A and B
    printf("Please enter the file path to matrix A (./examples/A_16x16): \n");
    fgets(userInputMatrixAFileName, 255, stdin);
    fflush(stdin);
    printf("\nPlease enter the file path to matrix B (./examples/B_16x16): \n");
    fgets(userInputMatrixBFileName, 255, stdin);
    fflush(stdin);

    char *matrixAFileName = strtok(userInputMatrixAFileName, "\n"),
         *matrixBFileName = strtok(userInputMatrixBFileName, "\n");

    printf("%s\n", matrixBFileName);
//    strcpy(matrixAFileName, "examples/A_16x16");
//    strcpy(matrixBFileName, "examples/B_16x16");

//    matrixAFileName = "./examples/A_16x16";
//    matrixBFileName = "./examples/B_16x16";
//    matrixAFileName = "./examples/A_64x64";
//    matrixBFileName = "./examples/B_64x64";

    // open matrix A file
    error = MPI_File_open(MPI_COMM_SELF, matrixAFileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    // exit on error
    if (error) {
        printf("Couldn't open matrix A file");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    // get matrix file size
    MPI_File_get_size(fh, &matrixFileSize);
    // calculate number of elements in the matrix
    totalElementCount = matrixFileSize / sizeof(double);

    // calculate matrix dimension
    matrixDimensions[activeMultCount - 1] = sqrt(totalElementCount);

    // define arrays to store matrices
    double matrixA[matrixDimensions[activeMultCount - 1]][matrixDimensions[activeMultCount - 1]],
           matrixB[matrixDimensions[activeMultCount - 1]][matrixDimensions[activeMultCount - 1]],
           *rcvMatrixA,
           *rcvMatrixB;

    // read matrix A
    error = MPI_File_read(fh, matrixA, matrixFileSize, MPI_DOUBLE, &status);
    if (error) MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_File_close(&fh);

    // open matrix B file
    error = MPI_File_open(MPI_COMM_SELF, matrixBFileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if (error) {
        printf("Couldn't open matrix B file");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    // read matrix B
    error = MPI_File_read(fh, matrixB, matrixFileSize, MPI_DOUBLE, &status);
    if (error) MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_File_close(&fh);

    // square root of amount of values per process
    int m = (int)sqrt(matrixDimensions[activeMultCount - 1]);

    double tmp;
    // determine number of workers
//     stop if it's not possible with at most 100 workers
    while(m > 0) {
//        printf("%d\n", matrixDimensions[activeMultCount - 1] * matrixDimensions[activeMultCount - 1]);
        tmp = sqrt(pow(matrixDimensions[activeMultCount - 1], 2) / pow(m, 2));
        if(tmp == (int)tmp && (int)sqrt(tmp) == sqrt(tmp)){
            workerCounts[activeMultCount - 1] = (int)tmp;
//            printf("%f\n", tmp);
//            printf("%d\n", (int)(workerCounts[activeMultCount - 1]));
            break;
        }
        m--;
    }

//    workerCounts[activeMultCount - 1] = 1;
    printf("%d workers will be spawned.\n", workerCounts[activeMultCount - 1]);

    // determine dimension of sub-matrices per process
    int subMatrixDimension = matrixDimensions[activeMultCount - 1]/sqrt(workerCounts[activeMultCount - 1]);

    // create new arguments array for worker processes
    char subMatrixDimensionStr[2],
         matrixDimension[2],
         **workerArgv;
    workerArgv = (char**)malloc(3 * sizeof(char*));

    sprintf(subMatrixDimensionStr, "%d", subMatrixDimension);    // convert subMatrixDimensionStr to string
    sprintf(matrixDimension, "%d", matrixDimensions[activeMultCount - 1]);          // convert matrixDimension to string

    workerArgv[0] = subMatrixDimensionStr;
    workerArgv[1] = matrixDimension;
    workerArgv[2] = NULL;

    // spawn workers
    MPI_Comm_spawn("./obj/worker.o", workerArgv, workerCounts[activeMultCount - 1], MPI_INFO_NULL, 0, MPI_COMM_SELF, &intercomms[activeMultCount - 1], MPI_ERRCODES_IGNORE);

    // create datatype for sub-matrices
    MPI_Datatype subMatrixDataType;
    int matrixSizes[2] = {matrixDimensions[activeMultCount - 1], matrixDimensions[activeMultCount - 1]},
        subMatrixSizes[2] = {subMatrixDimension, subMatrixDimension},
        starts[2] = {0, 0};

    MPI_Type_create_subarray(2, matrixSizes, subMatrixSizes, starts, MPI_ORDER_C, MPI_DOUBLE, &subMatrixDataType);

    MPI_Type_create_resized(subMatrixDataType, 0, subMatrixDimension * sizeof(double), &blockTypes[activeMultCount - 1]);

    // store new datatype in global blockTypes
    MPI_Type_commit(&blockTypes[activeMultCount - 1]);

    // distribute sub-matrices to workers per scatter
    int sendCountsPerProcess[workerCounts[activeMultCount - 1]],    // amount of units per process
        displacements[workerCounts[activeMultCount - 1]], // displacements per process
        tmpDisplacement = 0;

    for(i = 0; i < workerCounts[activeMultCount - 1]; i++){
        sendCountsPerProcess[i] = 1;
    }

    for(i = 0; i < sqrt(workerCounts[activeMultCount - 1]); i++){
        for(j = 0; j < sqrt(workerCounts[activeMultCount - 1]); j++){
            displacements[i * (int)sqrt(workerCounts[activeMultCount - 1]) + j] = tmpDisplacement;
//            printf("startMultiplications %d\n ", i * (int)sqrt(workerCounts[activeMultCount - 1]) + j);
            tmpDisplacement++;
        }
        tmpDisplacement += (subMatrixDimension - 1) * (int)sqrt(workerCounts[activeMultCount - 1]);
    }

//    for(i = 0; i < sqrt(workerCounts[activeMultCount - 1]); i++){
//        for(j = 0; j < sqrt(workerCounts[activeMultCount - 1]); j++){
//            printf("%d, ", displacements[i * (int)sqrt(workerCounts[activeMultCount - 1]) + j]);
//        }
//        printf("\n");
//    }

//    printf("\n%d\n", activeMultCount);

    // scatter matrix A
    MPI_Scatterv(
        matrixA,
        sendCountsPerProcess,
        displacements,
        blockTypes[activeMultCount - 1],
        rcvMatrixA,
        subMatrixDimension * subMatrixDimension,
//        pow(subMatrixDimension, 2),
        MPI_DOUBLE,
        MPI_ROOT,
        intercomms[activeMultCount - 1]
    );

    // scatter matrix B
    MPI_Scatterv(
        matrixB,
        sendCountsPerProcess,
        displacements,
        blockTypes[activeMultCount - 1],
        rcvMatrixB,
        subMatrixDimension * subMatrixDimension,
//        pow(subMatrixDimension, 2),
        MPI_DOUBLE,
        MPI_ROOT,
        intercomms[activeMultCount - 1]
    );

    // start listen for finishedMultiplications calculations
    MPI_Irecv(&finishedMultiplications, 1, MPI_INT, 0, 0, intercomms[activeMultCount - 1], &openRequests[activeMultCount - 1]);

//    free(matrixA);
}

void checkForRunningMultiplications(
    MPI_Comm intercomms[],
    MPI_Request openRequests[],
    MPI_Datatype blockTypes[],
    int *activeMultCount,
    int matricesDimensions[],
    int workerCounts[]
)
{
    MPI_Status status;
    int finishedFlag,
        currentMult,
        activeMultCountTmp = *activeMultCount,
        i;

    for(currentMult = activeMultCountTmp - 1; currentMult >= 0; currentMult--){
        MPI_Test(&openRequests[currentMult], &finishedFlag, &status);
        if(finishedFlag == 1){
            gatherResults(
                intercomms[currentMult],
                openRequests[currentMult],
                blockTypes[currentMult],
                &activeMultCountTmp,
                matricesDimensions[currentMult],
                workerCounts[currentMult],
                currentMult
            );
            for (i = currentMult; i < activeMultCountTmp; i++){
                intercomms[i] = intercomms[i + 1];
                openRequests[i] = openRequests[i + 1];
                blockTypes[i] = blockTypes[i + 1];
                workerCounts[i] = workerCounts[i + 1];
                matricesDimensions[i] = matricesDimensions[i + 1];
            }
            // decrease amount of running multiplications by finished multiplication
            *activeMultCount -= 1;
        }
    }
}


void waitForRunningMultiplications(
    MPI_Comm intercomms[],
    MPI_Request openRequests[],
    MPI_Datatype blockTypes[],
    int *activeMultCount,
    int matricesDimensions[],
    int workerCounts[]
)
{
    MPI_Status status;
    int finishedFlag,
        currentMult,
        activeMultCountTmp = *activeMultCount,
        i;

    for(currentMult = activeMultCountTmp - 1; currentMult >= 0; currentMult--){
        printf("Wait for multiplication %d to finish...\n", currentMult);
        MPI_Wait(&openRequests[currentMult], &status);
        gatherResults(
            intercomms[currentMult],
            openRequests[currentMult],
            blockTypes[currentMult],
            &activeMultCountTmp,
            matricesDimensions[currentMult],
            workerCounts[currentMult],
            currentMult
        );
        for(i = currentMult; i < activeMultCountTmp; i++){
            intercomms[i] = intercomms[i + 1];
            openRequests[i] = openRequests[i + 1];
            blockTypes[i] = blockTypes[i + 1];
            workerCounts[i] = workerCounts[i + 1];
            matricesDimensions[i] = matricesDimensions[i + 1];
        }
        // decrease amount of running multiplications by finished multiplication
        *activeMultCount--;
    }
}

void gatherResults(
    MPI_Comm intercomm,
    MPI_Request openRequest,
    MPI_Datatype blockType,
    int *activeMultCount,
    int matrixDimension,
    int workerCount,
    int currentMultId
){
    int sendCounts[workerCount],
        displacements[workerCount],
        subMatrixDimension = matrixDimension / sqrt(workerCount),
        tmpDisplacement = 0,
        i, j;
//    int displacements[4] = {0, 1, 16, 17};

    double matrixC[matrixDimension][matrixDimension],
           rcvMatrixC[subMatrixDimension][subMatrixDimension];

    char outputFileName[100];

    MPI_File fh;                    // file handle to open/save files

    for(i = 0; i < workerCount; i++){
        sendCounts[i] = 1;
    }

    for(i = 0; i < sqrt(workerCount); i++){
        for(j = 0; j < sqrt(workerCount); j++){
            displacements[i * (int)sqrt(workerCount) + j] = tmpDisplacement;
//            printf("gatherResults %d\n ", i * (int)sqrt(workerCount) + j);
            tmpDisplacement++;
        }
        tmpDisplacement += (subMatrixDimension - 1) * (int)sqrt(workerCount);
    }

//    printf("Displacements for %d subMatrixDimension: ", subMatrixDimension);
//    for(i = 0; i < 4; i++){
//        printf("%d, ", displacements[i]);
//    }
//    printf("\n");

    // collect multiplication result and save in matrix C
    MPI_Gatherv(
        rcvMatrixC,
        subMatrixDimension * subMatrixDimension,
        MPI_DOUBLE,
        matrixC,
        sendCounts,
        displacements,
        blockType,
        MPI_ROOT,
        intercomm
    );

    time_t ltime;       // calendar time
    ltime=time(NULL);   // get current cal time

    sprintf(outputFileName, "cannon_result_%d.dat", ltime);

    printf("Multiplication %d successfully finished!\n", currentMultId);
    printf("Results (complete for a 16x16 matrix):\n");
    for (i = 0; i < 16; i++){
        for (j = 0; j < 16; j++){
            printf("%.1f ", matrixC[i][j]);
            if(j == 15) printf("\n");
        }
    }
    printf("\n");

    // save result in file
    MPI_Status status;
    int error = MPI_File_open(MPI_COMM_SELF, outputFileName, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    if (error) MPI_Abort(MPI_COMM_WORLD, 1);
    error = MPI_File_write(fh, matrixC, pow(matrixDimension, 2), MPI_DOUBLE, &status);
    if (error) MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_File_close(&fh);

    printf("Result successfully saved to %s. \n", outputFileName);

    MPI_Type_free(&blockType);
}