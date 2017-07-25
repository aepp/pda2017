/***********************************************************************
Program:            Task2_master.c
Author:             Joshua Mendaza <joshua_mendaza@outlook.com>
matriclenumber:     4670232
Assignment:         5
Task:               2
Parameters:         int a, 1/0 decides if results are printed to the console
                    int secDelay, delays the worker processes for testing

Environment variables: no

Description:

The program listens to the input of the user.
The user can start a multiplication of 2 matrices A and B.
The filesnames of A and B are read via console input.
The master process spawns a group of workers, which multiply A and B
using the Cannon Algorithm.
There can be multiple multiplications at a time.
If the wokers finsihed, they send the result to the master,
who writes the result to a file.
**********************************************************************/

#include "mpi.h" 	    // import of the MPI definitions
#include <stdio.h> 	    // import of the definitions of the C IO library
#include <string.h>     // import of the definitions of the string operations
#include <unistd.h>	    // standard unix io library definitions and declarations
#include <errno.h>	    // system error numbers
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <sys/select.h>



/***********************************************************************
Name:           usage()
Parameters:     no
Return Value:   no

Function description:
The function is executed if -h is called as parameter
or not the right amount of parameters is passed to the program.
It prints an usage information for the program.
**********************************************************************/
void usage(int my_rank) {
    if (my_rank == 0){
        printf("\n"
        "usage: mpiexec -np 1 ./Task2_master printResult secDelay m\n"
        "\t-printResult \t1, if you want the result to be printed, 0 otherwise\n"
        "\t-secDelay \tamount of seconds all worker-groups are delayed (for testing)\n"
        "\nDescription\n"
        "\tThe program listens to the input of the user.\n"
        "\tThe user can start a multiplication of 2 matrices A and B.\n"
        "\tThe filesnames of A and B are read via console input.\n"
        "\tThe master process spawns a group of workers, which multiply A and B\n"
        "\tusing the Cannon Algorithm.\n"
        "\tThere can be multiple multiplications at a time.\n"
        "\tIf the wokers finsihed, they send the result to the master,\n"
        "\twho writes the result to a file.\n"
        "\nExample file names"
        "\n\tdata/A_16x16, data/B_16x16"
        "\n\tdata/A_64x64, data/B_64x64"
        "\n\tdata/A_256x256, data/B_256x256"
        "\n\nExample-call:\n"
        "\tmpiexec -np 1 ./Task2_master 1 30\n"
        "\n");
    }
}

/***********************************************************************
Name:           startMatrixMultiplication()
Parameters:     MPI_Comm intercomms[], array of intercomms, for the worker groups
                MPI_Request openRequests[], array of requests, to check if the workers are finish
                MPI_Datatype blocktypes[], blocktypes are used to scatter/gather the blocks of the matrices
                int matrixLengths[], size of matrices, for each multiplication
                int workerLengths[], amount of process (per row/column) for each multiplication
                int activeProcesses, amount of active mulitplications
Return Value:   no

Function description:
The function is executed if -h is called as parameter
or not the right amount of parameters is passed to the program.
It prints an usage information for the program.
**********************************************************************/
void startMatrixMultiplication(MPI_Comm intercomms[], MPI_Request openRequests[], MPI_Datatype blocktypes[], int matrixLengths[], int workerLengths[], int activeProcesses, int secDelay){

    int error;                          /// used to catch errors of functions
    char filenameMatrixARaw[255];       /// raw filename of matrix A
    char filenameMatrixBRaw[255];       /// raw filename of matrix B
    int matrixLength;                   /// length of rows and cols in the matrix
    int cellCount;                      /// number of cells in the matrix
    MPI_File fileMatrixA;               /// file of matrix A
    MPI_File fileMatrixB;               /// file of matrix B
    MPI_Status status;                  /// communication status


    /// get filename for matrix A from user
    printf("\nEnter filename for A:\n");
    fgets(filenameMatrixARaw, 255, stdin);

    fflush(stdin);

    /// get filename for matrix B from user
    printf("\nEnter filename for B:\n");
    fgets(filenameMatrixBRaw, 255, stdin);

    fflush(stdin);

    /// clean up the filenames (delete line breaks)
    char *filenameMatrixA = strtok(filenameMatrixARaw, "\n");
    char *filenameMatrixB = strtok(filenameMatrixBRaw, "\n");


    /// import Matrix A
    /// open file
    error = MPI_File_open(MPI_COMM_SELF, filenameMatrixA, MPI_MODE_RDONLY, MPI_INFO_NULL, &fileMatrixA);
    if (error) printf("Error 1");
    /// calculate size of matrix
    MPI_Offset MatrixASize;
    error = MPI_File_get_size(fileMatrixA, &MatrixASize);
    matrixLength = (int)sqrt((MatrixASize / sizeof(double)));
    cellCount = matrixLength * matrixLength;
    /// declare 2D arrays for matrix A
    double matrixA[matrixLength][matrixLength];
    double matrixB[matrixLength][matrixLength];
    /// read values from file
    error = MPI_File_read(fileMatrixA, matrixA, cellCount, MPI_DOUBLE, &status);
    if (error) printf("Error 2");
    /// close file
    MPI_File_close(&fileMatrixA);

    /// import Matrix B (same sizes are used)
    /// open file
    error = MPI_File_open(MPI_COMM_SELF, filenameMatrixB, MPI_MODE_RDONLY, MPI_INFO_NULL, &fileMatrixB);
    if (error) printf("Error 1");
    /// read values from file
    error = MPI_File_read(fileMatrixB, matrixB, cellCount, MPI_DOUBLE, &status);
    if (error) printf("Error 2");
    /// close file
    MPI_File_close(&fileMatrixB);


    /// calculate number of worker processes
    int a = (int)sqrt(sqrt(cellCount));     /// double square root of cells in the matrices
    int workerCount;                        /// total number of workers
    int workerLength;                       /// number of workers per row/column
    int tmp_workerCount;                    /// temporary workerCount
    int tmp_workerLength;                   /// temporary workerLength
    double cellsPerWorker;                  /// number of cells each worker gets
    double blockLengthPerWorker;            /// number of rows/cols each worker gets

    /// Basic idea:
    /// as default the double root of the cellCount is taken as workerCount
    /// in this case the blocksize of each process is the square root of cellCounts
    /// now we check if there is a smaller amount of workers which we can devide the cellCount
    /// there we iterate over the lenght of the process-matrix
    /// we take the first workerCount, which is smaller/equal to 100
    /// This limitation is done, to prevent a to small / to high amount of workers
    for (int i = a; i > 0; i--){

        tmp_workerLength = i;
        tmp_workerCount = i * i;

        cellsPerWorker = (double)cellCount / (double)tmp_workerCount;
        blockLengthPerWorker = (double)sqrt(cellsPerWorker);

        if (fmod(blockLengthPerWorker, 1) == 0.0){
            workerLength = tmp_workerLength;
            workerCount = tmp_workerCount;
        }
        if ((fmod(blockLengthPerWorker, 1) == 0.0) && (tmp_workerCount <= 100)){
            break;
        }
    }

    /// define MPI_Datatype for matrix blocks
    int rows = matrixLength;                /// row count of matrix
    int cols = matrixLength;                /// column count of matrix
    const int procRows = workerLength;      /// number of processes per column
    const int procCols = workerLength;      /// number of processes per row
    const int blockRows = rows / procRows;  /// number of rows in in one block (per process)
    const int blockCols = cols / procCols;  /// number of cols in in one block (per process)
    int blockSize = blockRows * blockCols;  /// number of cells per block
    MPI_Datatype blocktype2;                /// temporay blocktype

    /// the blocktype is stored in the array blocktypes, so we can reuse it later on
    MPI_Type_vector(blockRows, blockCols, cols, MPI_DOUBLE, &blocktype2);
    MPI_Type_create_resized( blocktype2, 0, sizeof(double), &blocktypes[activeProcesses-1]);
    MPI_Type_commit(&blocktypes[activeProcesses-1]);


    /// calculate the dicplacements for GatherV
    int disps[workerCount];
    int counts[workerCount];
    for (int ii=0; ii<procRows; ii++) {
        for (int jj=0; jj<procCols; jj++) {
            disps[ii*procCols+jj] = ii*cols*blockRows+jj*blockCols;
            counts [ii*procCols+jj] = 1;
        }
    }


    /// Spawn the worker process-group
    /// The intercommunicator is stored in the array, so we can use it later on
    MPI_Comm_spawn("./Task2_worker", MPI_ARGV_NULL, workerCount, MPI_INFO_NULL, 0, MPI_COMM_SELF, &intercomms[activeProcesses-1], MPI_ERRCODES_IGNORE);


    /// Broadcast relevant information to the workers
    /// BroadCast matrixLength, workerCount, workerLength, secDelay
    double blockA[blockRows][blockCols];    /// receive-matrix, not relevant for master-process
    double blockB[blockRows][blockCols];    /// receive-matrix, not relevant for master-process
    int metadata[4];                        /// send array, contains relevant metadata
    metadata[0] = matrixLength;
    metadata[1] = workerCount;
    metadata[2] = workerLength;
    metadata[3] = secDelay;
    /// Broadcast to all workers
    MPI_Bcast(metadata, 4, MPI_INT, MPI_ROOT, intercomms[activeProcesses-1]);

    /// Scatter Matrix A (by blocks) to all worker processes
    MPI_Scatterv(matrixA
                , counts
                , disps
                , blocktypes[activeProcesses-1]
                , blockA
                , blockRows*blockCols
                , MPI_DOUBLE
                , MPI_ROOT
                , intercomms[activeProcesses-1]);

    /// Scatter Matrix B (by blocks) to all worker processes
    MPI_Scatterv(matrixB
                , counts
                , disps
                , blocktypes[activeProcesses-1]
                , blockB
                , blockRows*blockCols
                , MPI_DOUBLE
                , MPI_ROOT
                , intercomms[activeProcesses-1]);


    /// store relevant information in the arrays, so we can use it later on
    matrixLengths[activeProcesses-1] = matrixLength;
    workerLengths[activeProcesses-1] = workerLength;

    /// Set a nonBlocking receive to one worker process
    /// This acts as action listener, we can test, if the mulitplication finished
    /// The request is stored in the request array, so we can use it later on
    int tmp;
    MPI_Irecv(&tmp, 1, MPI_INT, 0, 1, intercomms[activeProcesses-1], &openRequests[activeProcesses-1]);

    /// Print info to user
    printf("\n\nMatrix-Multiplication started succesfully"
    "\nworkerCount: %d, workerLength: %d", workerCount, workerLength);

}



/***********************************************************************
Name:           testRequest()
Parameters:     MPI_Request request
Return Value:   int finished

Function description:
The function tests a the request, of a non-blocking communciation.
It returns the status of the request.
**********************************************************************/
int testRequest(MPI_Request request){
    int finished;
    MPI_Status status;
    MPI_Test(&request, &finished, &status);
    return finished;
}

/***********************************************************************
Name:           rem***FromArray()
Parameters:     array[], of different datatype
                int ***Count, number of entries in the array
                int ***ToDel, index of entries, which should be deleted
Return Value:   no

Function description:
The functions delete an entrie from an array.
Therefore the elements with higher index are shifted to the left.
**********************************************************************/
void remReqFromArray(MPI_Request requests[], int reqCount, int reqToDel){
    for (int i = reqToDel; i < reqCount-1; i++){
        requests[i] = requests[i+1];
    }
}


void remComFromArray(MPI_Comm comms[], int comCount, int comToDel){
    for (int i = comToDel; i < comCount-1; i++){
        comms[i] = comms[i+1];
    }
}


void remDatFromArray(MPI_Datatype dats[], int datCount, int datToDel){
    for (int i = datToDel; i < datCount-1; i++){
        dats[i] = dats[i+1];
    }
}


void remIntFromArray(int ints[], int intCount, int intToDel){
    for (int i = intToDel; i < intCount-1; i++){
        ints[i] = ints[i+1];
    }
}



/***********************************************************************
Name:           checkStatus()
Parameters:     MPI_Comm intercomms[], array of intercomms, for the worker groups
                MPI_Request openRequests[], array of requests, to check if the workers are finish
                MPI_Datatype blocktypes[], blocktypes are used to scatter/gather the blocks of the matrices
                int matrixLengths[], size of matrices, for each multiplication
                int workerLengths[], amount of process (per row/column) for each multiplication
                int *activeProcesses, amount of active mulitplications
                int *finishedProcesses, amount of multiplications which are finished
                int finishAll, flag if 1, then we wait for all open mulitplications to finish
                int printResult, flag if 1, then we print the result
Return Value:   no

Function description:
The function checks the status for each worker-group.
If the workers are finish, the master process receives the result
and writes it to a file.
**********************************************************************/
void checkStatus(MPI_Comm intercomms[], MPI_Request openRequests[], MPI_Datatype blocktypes[], int matrixLengths[], int workerLengths[], int *activeProcesses, int *finishedProcesses, int finishAll, int printResult){
    MPI_Status status;                  /// communication status
    int actProBuf = *activeProcesses;   /// tmp copy of the number of active processes

    /// iterate over all active multiplication
    /// we iterate from highest index to lowest, because we delete from arrays
    for (int activeProcess = actProBuf-1; activeProcess >= 0; activeProcess--){

        /// get metadata of the multiplication, by accessing the arrays
        int matrixLength = matrixLengths[activeProcess];
        int workerLength = workerLengths[activeProcess];
        int workerCount = workerLength * workerLength;

        int finished;                   /// status of the worker group


        /// if the program wasnt quited by the user, we just check if the multiplcation finished
        /// if it isn't finish, we skip it
        if (finishAll == 0){
            /// we use the non-blocking receive, to check if the worker group finished
            finished = testRequest(openRequests[activeProcess]);
        // if programm was quited, than we wait for the workergroup to finish
        }else {
            printf("\n\nWait for Run: %d", activeProcess);
            MPI_Wait(&openRequests[activeProcess], &status);
            finished = 1;
        }

        /// if the worker group finished
        if (finished != 0){
            /// define displacements and counts for gather v
            int rows = matrixLength;                /// row count of matrix
            int cols = matrixLength;                /// column count of matrix
            const int procRows = workerLength;      /// number of processes per column
            const int procCols = workerLength;      /// number of processes per row
            const int blockRows = rows / procRows;  /// number of rows in in one block (per process)
            const int blockCols = cols / procCols;  /// number of cols in in one block (per process)
            int blockSize = blockRows * blockCols;  /// number of cells per block

            int disps[workerCount];
            int counts[workerCount];
            for (int ii=0; ii<procRows; ii++) {
                for (int jj=0; jj<procCols; jj++) {
                    disps[ii*procCols+jj] = ii*cols*blockRows+jj*blockCols;
                    counts [ii*procCols+jj] = 1;
                }
            }

            double resultC[matrixLength][matrixLength];     /// 2D array for the result
            double send[blockRows][blockCols];              /// receive array, used by the workers
            MPI_Gatherv(send, blockSize, MPI_DOUBLE, resultC, counts, disps, blocktypes[activeProcess], MPI_ROOT, intercomms[activeProcess]);


            /// print the result (if the user wants it to be printed)
            if (printResult == 1){
                for (int i = 0; i < matrixLength; i++){
                    printf("\n");
                    for (int j = 0; j < matrixLength; j++){
                        printf("%lf, ", resultC[i][j]);
                    }
                }
            }


            /// Write the result to the file
            /// The filename is a consequtive number, which is declared by finishedProcesses
            *finishedProcesses += 1;
            char resultFileName[100];
            int aInt = *finishedProcesses;
            char str[15];
            sprintf(str, "%d", aInt);
            memcpy(resultFileName, "result_", 7 * sizeof(char));
            memcpy(resultFileName+7, str, 15 * sizeof(char));

            /// open result file
            MPI_File resultFile;
            MPI_File_open(MPI_COMM_SELF, resultFileName, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &resultFile);

            /// write values to the result file
            MPI_File_write_ordered(resultFile, resultC, (matrixLength * matrixLength), MPI_DOUBLE, &status);

            /// close result file
            MPI_File_close(&resultFile);

            /// clean up arrays
            /// This part is important!!!
            /// if we don't clean up, than a worker-group, which is already finish, would be checked
            *activeProcesses -= 1;
            remReqFromArray(openRequests, actProBuf, activeProcess);
            remComFromArray(intercomms, actProBuf, activeProcess);
            remDatFromArray(blocktypes, actProBuf, activeProcess);
            remIntFromArray(matrixLengths, actProBuf, activeProcess);
            remIntFromArray(workerLengths, actProBuf, activeProcess);
            printf("\nRun %d finished and written to: %s", activeProcess, resultFileName);

        } else {
            printf("\nRun %d not yet finish!", activeProcess);
        }
    }
}


/***********************************************************************
Name:           main()
Parameters:     argc number of parameters passed to the program
                argv array of parameters passed to the program
Return Value:   int

Function description:
The main function of the program.
**********************************************************************/
int main(int argc, char* argv[ ])
{

    /// Declaration of variables
	int process_count;                  // amount of running processes
	int namelen;                        // length of name
	int my_rank;                        // rank of the process in MPI_COMM_WORLD
	char *c, proc_name[MPI_MAX_PROCESSOR_NAME+1]; 	// hostname
    double start_time;                  // timestamp after MPI_Init()
    int error;                          // return value of various function-calls

    /// Initializing of MPI-Interface
	MPI_Init(&argc, &argv);
    /// meassure start-time
    start_time = MPI_Wtime();
    /// Get amount of running processes
	MPI_Comm_size(MPI_COMM_WORLD, &process_count);
    /// Get your rank
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    /// Initialize string with NULL-characters
	memset (proc_name, 0, MPI_MAX_PROCESSOR_NAME+1 );
    /// Finding out own computer name
	MPI_Get_processor_name(proc_name, &namelen);
    /// Separate the first part of hostname
	if ( (c=strchr(proc_name,'.'))  !=  NULL) *c = '\0';

    /// Check if right amount of parameters exist
    if((argc != 3) || (strcmp(argv[1], "-h") == 0)){
        usage(my_rank);
        MPI_Finalize();
        return -1;
    }

    int printResult = atoi(argv[1]);        // 1 = results are printed, 0 = results are not printed
    int secDelay = atoi(argv[2]);           // delays the workers for x seconds
    int maxParaMulti = 10;                   // Max number of mulitplications at the same time


    int activeProcesses = 0;                /// number of active mulitplications
    int totalProcesses = 0;                 /// number of total mulitplications
    int finishedProcesses = 0;              /// number of finished mulitplications
    MPI_Request openRequests[maxParaMulti]; /// array of requests, to check if the workers are finish
    MPI_Comm intercomms[maxParaMulti];      /// array of intercomms, for the worker groups
    MPI_Datatype blocktypes[maxParaMulti];  /// blocktypes are used to scatter/gather the blocks of the matrices
    int matrixLengths[maxParaMulti];        /// size of matrices, for each multiplication
    int workerLengths[maxParaMulti];        /// amount of process (per row/column) for each multiplication


    /// Set up loop for the user interaction

    int stop = 0;           /// flag if user interaction should be stopped
    int fd = 0;             /// file discriptor for standard-input
    char buf[11];           /// input buffer
    int ret;                /// return of the read function
    int sret;               /// return of the select function
    fd_set readfs;          /// array of file-discriptors for the select function
    struct timeval timeout; /// timeout for the select function

    /// loop until the program should be stopped
    while(stop == 0){

        /// add standard input to the select function
        FD_ZERO(&readfs);
        FD_SET(fd, &readfs);

        /// set the input time to 5 seconds
        /// we wait 5 seconds for the user to input something
        timeout.tv_sec = 5;
        timeout.tv_usec = 0;

        /// print info to user
        printf("\n\nWhat do u wanna do? (s)tart, (q)uit\n");
        /// start select functions
        sret = select(8, &readfs, NULL, NULL, &timeout);

        /// handle return values of the select function

        /// if the user made an input
        if (sret != 0){
            /// read the input
            memset((void *) buf, 0, 11);
            ret = read(fd, (void *) buf, 10);

            /// if the input startet with an 's'
            if (buf[0] == 's'){

                ///Check if we reached maxParaMulti
                if (activeProcesses == maxParaMulti){
                    printf("\nSorry, you have already %d multiplcations running."
                    "\nPlease wait for one to finish, to start another one.\n", maxParaMulti);
                } else{
                    /// increase the number of active multiplications
                    activeProcesses += 1;
                    totalProcesses += 1;
                    /// Start the matrix-multiplication
                    startMatrixMultiplication(intercomms, openRequests, blocktypes, matrixLengths, workerLengths, activeProcesses, secDelay);
                }

            /// if the input started with 'q' than stop the loop
            } else if(buf[0] == 'q'){
                stop = 1;

            /// in case (wrong input), do nothing
            } else{
                stop = 0;
            }

        /// if the user made no input (after 5 seconds)
        } else{

            /// if running multiplcations, than check their status
            if(activeProcesses > 0){
                printf("\nCheck running multiplications ...\n");
                checkStatus(intercomms, openRequests, blocktypes, matrixLengths, workerLengths, &activeProcesses, &finishedProcesses, 0, printResult);
            } else{
                printf("\nNo running mulitplications\n");
            }
        }

    }

    /// after the loop was quitted, we wait for all worker-groups to finish
    printf("\nNow, we wait for all processes to finish ...\n");
    checkStatus(intercomms, openRequests, blocktypes, matrixLengths, workerLengths, &activeProcesses, &finishedProcesses, 1, printResult);


    MPI_Finalize();
    return 0;
}
