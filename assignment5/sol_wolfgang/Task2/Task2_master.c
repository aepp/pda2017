/************************************************************************/
/* Program: Task2_master.c                                              */
/* Author: Wolfgang Stammer <wolfgang.stammer@stud.uni-frankfurt.de>    */
/* matriclenumber: 5962350                                              */
/* Assignment : 5                                                       */
/* Task: 2                                                              */
/* Parameters:  -h: view help                                           */
/*                                                                      */
/* Environment variables: no                                            */
/*                                                                      */
/* Description:                                                         */
/* This Program performs Task 2 of Assignment 5. It begins a master     */
/* process which starts the select function in order to listen to the   */
/* standard input from command line. The user is then asked whether to  */
/* start a multiplication of two nxn matrices or exit the program. If   */
/* start is entered a new multiplication is started and solved using    */
/* the Cannon algorithm on several worker processes. The matrix files   */
/* and the number of worker processes is given by the user during       */
/* runtime. Once the multiplication is started on the worker processes  */
/* the user is again asked whether to start a multiplication or exit    */
/* the program. If the user enters nothing the question is repeated     */
/* every 5 seconds and internally it is checked if an active            */
/* multiplication is finished and can be written to a file. A new       */
/* multiplication can be started when ever the user is asked again.     */
/* However there is a maximum of possible running multiplications, such */
/* that the user must first wait for an active run to finish, before    */
/* starting a new one. This ensures that the cluster is not overloaded. */
/* If the user request to exit, all runs are tested if they have        */
/* finished and if not the program waits for all runs to finish, before */
/* receiving the result matrix C from the worker processes and saving   */
/* them to a new file.                                                  */
/*                                                                      */
/************************************************************************/

#include "mpi.h" 	    // import of the MPI definitions
#include <stdio.h> 	    // import of the definitions of the C IO library
#include <time.h>       // import rand seed generator, elapsed time measurement
#include <stdlib.h>
#include <string.h>     // import of the definitions of the string operations
#include <unistd.h>	    // standard unix io library definitions and declarations
#include <errno.h>	    // system error numbers
#include <math.h>
#include <linux/limits.h>   // for max path length
#include <sys/select.h>

void rem_Elem_from_Array(MPI_Comm [], MPI_Request [], MPI_Datatype [], int [],
                            int [], int [], int [], int, int);
void checkRequest(MPI_Request, int *);
void checkRuns(MPI_Comm [], MPI_Request [], MPI_Datatype [], int [], int [],
                int [], int [], int *, int);
void startMult(MPI_Comm [], MPI_Request [], MPI_Datatype [], int [], int [],
                int [], int, int);
void init_ones(int, int[]);
void get_displacements(int, int, int []);
int check_help(int, char**, int, int);
void help_Info();

int main(int argc, char* argv[ ]){

    const int root = 0;
    int     my_rank, 					// rank of the process
            numprocs,                   // number of processes
            fd = 0,                     // file discriptor
            ret,                        // return of read function
            sel_ret,                    // select return value
            loop_stop = 0,
            totalRuns = 0,           // indicator if any mult. was started
            activeRuns = 0,            // number of active multiplicaitons
            maxNumMult = 5,             // maximal number of mult that can be
                                        // started, to prevent overload
            nps[maxNumMult],            // array of number of worker processes
            dim_Mats[maxNumMult],       // array storing the original matrix
                                        // dimensions for all runs
            dim_SubMats[maxNumMult],    // array storing the submatrix
                                        // dimensions for all runs
            globalRunIDs[maxNumMult];   // array storing the global run index,
                                        // used for file writing
    char    buf[11];                    // stdin buffer
    fd_set  readfds;                    // read file desctiptor socket for select()
    struct timeval timeout;             // wait time of select()
    MPI_Comm intercomms[maxNumMult];    // array of intercommunicators
    MPI_Request openRequests[maxNumMult];   // array of requests
    MPI_Datatype blocktypes[maxNumMult];    // array of sub matrix datatypes

    /* check for help request and if necessary quit program */
    if (check_help(argc, argv, my_rank, root) == -1){
        return -1;
    }

    /* initializing of MPI-Interface */
    MPI_Init(&argc, &argv);
    /* get rank IDs */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    while (loop_stop == 0){

        // add stdin to select function
        FD_ZERO(&readfds);
        FD_SET(fd, &readfds);

        // set time that program initially waits for input to 5 seconds
        timeout.tv_sec = 5;
        timeout.tv_usec = 0;

        // ask user for input
        printf("\nWhat would you like to do? Please enter 'start' or 'quit':\n");
        // start select()
        sel_ret = select(1, &readfds, NULL, NULL, &timeout);

        // if input was made
        if(sel_ret !=0){

            // read input
            memset((void *) buf, 0, 11);
            ret = read(fd, (void *) buf, 10);
            // remove \0 from string for strcmp
            buf[strlen(buf)-1] = 0;

            // differentiate between start and quit entries
            if (strcmp(buf, "start") == 0){

                printf("You have entered: %s\n", buf);

                // if maximum number of multiplications reached
                if (activeRuns == maxNumMult){
                    printf("\nSorry you have already reached the maximum "\
                            "number of parallel \nmultiplications that can run "\
                            "simultaneously (%d).\nPlease wait for one to "\
                            "finish.\n", maxNumMult);
                }
                // otherwise start new multiplication
                else {
                    printf("\nStarting new multiplication!\n");
                    /* increase both count variables */
                    activeRuns += 1;
                    totalRuns += 1;
                    /* add new run to array of global IDS */
                    globalRunIDs[activeRuns-1] = totalRuns;
                    /* start multiplication run */
                    startMult(intercomms, openRequests, blocktypes, nps,
                                dim_Mats, dim_SubMats, activeRuns, root);
                }
            }
            // otherwise quit program
            else if ((strcmp(buf, "quit") == 0)){
                printf("You have entered: %s\n", buf);
                printf("The program will end now.\n");
                loop_stop = 1;
            }
        }
        // if no input was made
        else {
            // check how many calculations are still running
            if (activeRuns > 0){
                printf("\n%d multiplication(s) still running...\n", activeRuns);
                checkRuns(intercomms, openRequests, blocktypes, nps, dim_Mats,
                            dim_SubMats, globalRunIDs, &activeRuns, 0);
            }
            else {
                printf("\nThere is no running multiplication...\n");
            }
        }
    }
    if (activeRuns != 0){
        printf("\nWaiting for all processes to finish...\n");
        checkRuns(intercomms, openRequests, blocktypes, nps, dim_Mats,
                    dim_SubMats, globalRunIDs, &activeRuns, 1);
        printf("\n");
    }
    else printf("\n");

    MPI_Finalize();		            // finalizing MPI interface

    return 0;
}

/*******************************************************************************
Name:       rem_Elem_from_Array
Parameters: MPI_Comm intercomms[]: some array
            MPI_Request openRequests[]: some array
            MPI_Datatype blocktypes[]: some array
            int nps[]: some array
            int dim_SubMats[]: some array
            int dim_Mats[]: some array
            int globalRunIDs[]: some array
            int IDX_toRem: index of element to be removed form arrays
            int IDX_lastElem: how many elements are stored in each array
Return:

Description:
This function receives various arrays of different types and removes the
element from each array at index IDX_toRem. This is done by shifting all
elements after the element at index IDX_toRem one to the left.
*******************************************************************************/
void rem_Elem_from_Array(MPI_Comm intercomms[], MPI_Request openRequests[],
                        MPI_Datatype blocktypes[], int nps[], int dim_SubMats[],
                        int dim_Mats[], int globalRunIDs[], int IDX_toRem,
                        int IDX_lastElem){
    int i;
    for (i = IDX_toRem; i < IDX_lastElem; i++){
        intercomms[i] = intercomms[i+1];
        openRequests[i] = openRequests[i+1];
        blocktypes[i] = blocktypes[i+1];
        nps[i] = nps[i+1];
        dim_SubMats[i] = dim_SubMats[i+1];
        dim_Mats[i] = dim_Mats[i+1];
        globalRunIDs[i] = globalRunIDs[i+1];
    }

}

/*******************************************************************************
Name:       checkRequest
Parameters: MPI_Request request: request of current run
            int *finished: integer pointer, stores result of MPI_Test()
Return:

Description:
This function receives an MPI_Request, request, of the current run and tests
whether the worker processes have sent to the master process that the
calculations are finished. The test is done using MPI_Test() and the result is
stored in the pointer finished.
*******************************************************************************/
void checkRequest(MPI_Request request, int *finished){
    MPI_Status status;
    MPI_Test(&request, finished, &status);
}

/*******************************************************************************
Name:       checkRuns
Parameters: MPI_Comm intercomms[]: array of intercommunicators of active runs
            MPI_Request openRequests[]: array of finish request of active runs
            MPI_Datatype blocktypes[]: array of blocktypes of active runs
            int nps[]: array of number of worker processes of active runs
            int dim_Mats[]: array of original matrix dimension of active runs
            int dim_SubMats[]: array of submatrix dimensions of active runs
            int globalRunIDs[]: array of global IDs of active runs
            int *activeRuns: integer of number active runs
            int finishAll: integer indicator whether all runs should be finished
Return:

Description:
This function recieves arrays of information of all active multiplication runs,
namely intercomms, openRequests, blocktypes, nps, dim_Mats, dim_SubMats and
globalRunIDs. These contain the intercommunicators, finish requests, number of
processes, original matrix dimensions, submatrix dimension of each worker
process, and global run ID, of each multiplication run, in the listed order.
The function goes through all active runs from last to first and checks if the
workers of the run have finished calculations, using checkRequest(). If all
runs should be finished, indicated by finishAll, and the current run is not
finished, the function waits as long as it takes for the run to be finished.
Once the run is finished MPI_Gatherv() is used to collect the data from the
workers to the master process, using specific displacements and counts. A new
file is opened and the result matrix is written to the file. Also the first
10 columns of the first 10 rows are printed. The information from the current
run in all relevant arrays are deleted.
*******************************************************************************/
void checkRuns(MPI_Comm intercomms[], MPI_Request openRequests[],
                        MPI_Datatype blocktypes[], int nps[], int dim_Mats[],
                        int dim_SubMats[], int globalRunIDs[], int *activeRuns,
                        int finishAll){

    int     activeRunsIDX,                  // index of current run
            tmp_activeRuns = *activeRuns,   // number of running mult. before
                                            // finishing any
            err = 0,
            i, j,
            finished,                   // indicator if run is finished
            dim_Mat,            // dimension of original matrix of current run
            dim_SubMat;         // dimension of submatrix of current run
    char    fileNameOut[22];            // file path for results
    double  *matrixC,                   // final result array pointer
            *rcv_matrixC;               // pointer, not significant for root
    MPI_File fh;                        // file handle
    MPI_Status status;

    /* go through all running multiplications and check if finished */
    for (activeRunsIDX = tmp_activeRuns-1; activeRunsIDX >= 0; activeRunsIDX--){

        /* if not all calculations should be finished just test if current */
        /* calculation is finished */
        /* strangely 39 indicates complete, 59 indicates not complete*/
        checkRequest(openRequests[activeRunsIDX], &finished);

        /* otherwise wait for the calculation to finish */
        if((finished != 1) && (finishAll == 1)){
            MPI_Wait(&openRequests[activeRunsIDX], &status);
            finished = 1;
        }

        /* all worker processes of calculation finished */
        /* gather the results and save to file */
        if (finished == 1) {
            /* get dimensions for this calculation */
            dim_Mat = dim_Mats[activeRunsIDX];
            dim_SubMat = dim_SubMats[activeRunsIDX];

            /* set memory for result array */
            matrixC = malloc(dim_Mat * dim_Mat * sizeof(double));

            // how many pieces of data everyone has, in units of blocks
            int sndcounts[nps[activeRunsIDX]];
            init_ones(nps[activeRunsIDX], sndcounts);
            // the starting point of everyone's data in the global array, in block extents
            int snddispls[nps[activeRunsIDX]];
            get_displacements(nps[activeRunsIDX], dim_SubMat, snddispls);

            /* gather all computed results from each worker process */
            MPI_Gatherv(rcv_matrixC, dim_SubMat*dim_SubMat,  MPI_DOUBLE,
                        matrixC, sndcounts, snddispls, blocktypes[activeRunsIDX],
                        MPI_ROOT, intercomms[activeRunsIDX]);

            /* create name for output file */
            sprintf(fileNameOut, "Result_MatrixC_%d.dat", globalRunIDs[activeRunsIDX]);

            /* write result matrix to output file */
            err = MPI_File_open(MPI_COMM_SELF, fileNameOut, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
            if (err) MPI_Abort(MPI_COMM_WORLD, 1);  // exit on error
            err = MPI_File_write(fh, matrixC, dim_Mat*dim_Mat, MPI_DOUBLE, &status);
            if (err) MPI_Abort(MPI_COMM_WORLD, 1);  // exit on error
            MPI_File_close(&fh);

            printf("\nThe resulting matrix C of run %d has been saved to the "\
                    "file %s. \n", globalRunIDs[activeRunsIDX] , fileNameOut);

            printf("\nHere are the first 10 columns of the first 10 rows of the resulting\n"
                    "matrix of run %d: \n \n", globalRunIDs[activeRunsIDX]);
            for (i = 0; i < 16; i++){
                for (j = 0; j < 16; j++){
                    printf("%lf ", matrixC[i*dim_Mat+ j]);
                    if(j == 15) printf("\n");
                }
            }

            /* free datatype */
            MPI_Type_free(&blocktypes[activeRunsIDX]);

            /* free memory */
            free(matrixC);

            /* remove current entry in all relevant arrays */
            rem_Elem_from_Array(intercomms, openRequests, blocktypes, nps,
                                    dim_SubMats, dim_Mats, globalRunIDs,
                                    activeRunsIDX, tmp_activeRuns);

            /* reduce number of running calculations */
            *activeRuns -= 1;

        }
        // otherwise run isn't finished
        else{
            printf("\nRun %d has not yet finished calculations.\n", globalRunIDs[activeRunsIDX]);
        }
    }

}

/*******************************************************************************
Name:       startMult
Parameters: MPI_Comm intercomms[]: array of intercommunicators of active runs
            MPI_Request openRequests[]: array of finish request of active runs
            MPI_Datatype blocktypes[]: array of blocktypes of active runs
            int nps[]: array of number of worker processes of active runs
            int dim_Mats[]: array of original matrix dimension of active runs
            int dim_SubMats[]: array of submatrix dimensions of active runs
            int activeRuns: integer of number active runs
            int root: integer of root process ID
Return:

Description:
This function starts the Cannon multiplication algorithm on several worker
processes. It receives the file paths for Matrix A and Matrix B and the number
of worker processes. If the number of worker processes does not fit to the size
of the matrices the user is prompted again for the number of worker processes
(for details see code). The matrix files are then opened and read by the master
process. The master process then spawns the worker processes and then sends a
submatrix of both arrays to each worker process using MPI_Scatterv() using
a displacements array and count array. Finally a nonblocking communication is
started to listen for when all worker processes have finished calculations. All
new informations for the newly started multiplication run are stored in the
relevant arrays (intercomms, openRequests, blocktypes, nps, dim_Mats,
dim_SubMats).
*******************************************************************************/
void startMult(MPI_Comm intercomms[], MPI_Request openRequests[],
                        MPI_Datatype blocktypes[], int nps[], int dim_Mats[],
                        int dim_SubMats[], int activeRuns,
                        int root){

    int     dim_Mat,                    // dimension of matrices
            dim_SubMat,                 // dimension of submatrices of workers
            np = 0,                     // number of worker processes
            err = 0,
            i,j,
            finished;
    double  *matrixA,                   // matrix A pointer array
            *rcv_matrixA,               // not significant for root
            *matrixB,                   // matrix B pointer array
            *rcv_matrixB;               // not significant for root
    char    *fileNameA,        // file path for Matrix A
            *fileNameB;        // file path for Matrix B
    MPI_File fh;                        // file handle
    MPI_Status status;
    MPI_Offset fileSize;

    /* set memory for file path variables */
//    memset(fileNameA, 0, PATH_MAX);
//    memset(fileNameB, 0, PATH_MAX);

    /* receive file paths for matrices A and B */
    // disable output buffer
//    setbuf(stdout, NULL);
//    printf("\nPlease enter the file path to matrix A, e.g. '5/A_16x16' : \n \n");
//    scanf("%s", fileNameA);
//    printf("\nThe file path to matrix A is: %s\n", fileNameA);
//
//    printf("\nPlease enter the file path to matrix B, e.g. '5/B_16x16' : \n \n");
//    scanf("%s", fileNameB);
//    printf("\nThe file path to matrix B is: %s\n", fileNameB);
    fileNameA = "5/A_16x16";
    fileNameB = "5/B_16x16";

    /* Open MatrixA file */
    err = MPI_File_open(MPI_COMM_SELF, fileNameA, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if (err) MPI_Abort(MPI_COMM_WORLD, 1);  // exit on error
    /* get size of matrix */
    MPI_File_get_size(fh, &fileSize);
    fileSize = fileSize/sizeof(double);     // number of pixels altogether
    // dimension of matrices in either direction
    dim_Mat = sqrt(fileSize);
    dim_Mats[activeRuns-1] = dim_Mat;

    /* allocate matrixA and matrixB memory */
    matrixA = malloc(dim_Mat * dim_Mat * sizeof(double));
    matrixB = malloc(dim_Mat * dim_Mat * sizeof(double));

    /* read matrixA from file */
    err = MPI_File_read(fh, matrixA, fileSize, MPI_DOUBLE, &status);
    if (err) MPI_Abort(MPI_COMM_WORLD, 1);  // exit on error
    MPI_File_close(&fh);                    // close matrix A file

    /* Open MatrixB file */
    err = MPI_File_open(MPI_COMM_SELF, fileNameB, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if (err) MPI_Abort(MPI_COMM_WORLD, 1);  // exit on error
    /* read matrixB from file */
    err = MPI_File_read(fh, matrixB, fileSize, MPI_DOUBLE, &status);
    if (err) MPI_Abort(MPI_COMM_WORLD, 1);  // exit on error
    MPI_File_close(&fh);                    // close matrix A file

    /* receive number of worker processes from user, make sure proper number */
    /* i.e. n/sqrt(p) is whole number */
    do{
//        printf("\nPlease enter the number of worker processes:\n");
//        scanf("%d", &np);
        np = 4;
        if ((int) (dim_Mat/sqrt(np)) != dim_Mat/sqrt(np)){
            printf("\nERROR: The number of rows must be equally dividable by the number of processes! \n"
                    "The number of rows for the specified matrices is %d. Please restart the \n"
                    "program with a different number of processes. \n ", dim_Mat);
        }
    }while((int) (dim_Mat/sqrt(np)) != dim_Mat/sqrt(np));
    printf("\nThe number of processes is: %d\n", np);
    // store number of workers
    nps[activeRuns-1] = np;

    /* compute dimension of submatrices */
    dim_SubMat = dim_Mat/sqrt(np);
    dim_SubMats[activeRuns-1] = dim_SubMat;

    /* create new argv array with relevant variables for worker processes */
    char str1[2];
    char str2[2];
    sprintf(str1, "%d", dim_SubMat);        // convert dim_SubMat to string
    sprintf(str2, "%d", dim_Mat);           // convert dim_Mat to string
    char **new_argv;
    new_argv = (char**)malloc(3*sizeof(char*));
    new_argv[0] = str1;
    new_argv[1] = str2;
    new_argv[2] = NULL;
    /* spawn worker processes */
    MPI_Comm_spawn("./Task2_workers", new_argv, np, MPI_INFO_NULL, root, MPI_COMM_SELF,\
                    &intercomms[activeRuns-1], MPI_ERRCODES_IGNORE);

    /* create submatrix datatype */
    MPI_Datatype submat_type;
    int sizes[2] = {dim_Mat, dim_Mat},          /* size of entire matrix */
        subsizes[2] = {dim_SubMat, dim_SubMat}, /* size of sub-matrix */
        starts[2] = {0,0};
    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &submat_type);
    /* change the extent of the type */
    MPI_Type_create_resized(submat_type, 0, dim_SubMat*sizeof(double), &blocktypes[activeRuns-1]);
    // save in blocktype array
    MPI_Type_commit(&blocktypes[activeRuns-1]);

    /* scatter submatrices from root to all worker processes */
    // how many pieces of data everyone has, in units of blocks
    int sndcounts[np];
    init_ones(np, sndcounts);
    // the starting point of everyone's data in the global array, in block extents
    int snddispls[np];
    get_displacements(np, dim_SubMat, snddispls);

//    for(i = 0; i < sqrt(np); i++){
//        for(j = 0; j < sqrt(np); j++){
//            printf("%d, ", snddispls[i * (int)sqrt(np) + j]);
//        }
//        printf("\n");
//    }
//    printf("\n%d\n", activeRuns);

    // scatter matrix A
    MPI_Scatterv(
        matrixA,
        sndcounts,
        snddispls, /* proc i gets counts[i] types from displs[i] */
        blocktypes[activeRuns-1],
        rcv_matrixA,
        dim_SubMat*dim_SubMat,
        MPI_DOUBLE,
        MPI_ROOT,
        intercomms[activeRuns-1]
    );

    // scatter matrix B
    MPI_Scatterv(matrixB, sndcounts, snddispls, /* proc i gets counts[i] types from displs[i] */
                blocktypes[activeRuns-1], rcv_matrixB, dim_SubMat*dim_SubMat, MPI_DOUBLE,
                MPI_ROOT, intercomms[activeRuns-1]);

    // initiate non-blocking communication to process 0 of worker processes
    // to listen for finished calculations
    MPI_Irecv(&finished, 1, MPI_INT, 0, 0, intercomms[activeRuns-1], &openRequests[activeRuns-1]);

    /* free memory */
    free(new_argv);
    free(matrixA);
    free(matrixB);
}

/*******************************************************************************
Name:       init_ones
Parameters: int size: size of array
            int array[size]: array of integers
Return:

Description:
This function receives an integer array and writes it with ones.
*******************************************************************************/
void init_ones(int size, int array[size]){
    int i;
    for (i = 0; i < size; i++){
        array[i] = 1;
    }
}

/*******************************************************************************
Name:       get_displacements
Parameters: int np: number of processes
            int dim: dimension of blocks
            int array[np]: array with displacements stored in it
Return:

Description:
This function computes the displacements of np blocks with dimension dim by
block extent and saves the values in array.
*******************************************************************************/
void get_displacements(int np, int dim, int array[np]){
    int disp = 0;
    for (int i=0; i<sqrt(np); i++) {
        for (int j=0; j<sqrt(np); j++) {
            array[i*(int)sqrt(np)+j] = disp;
            printf("%d, ", i*(int)sqrt(np)+j);
            disp += 1;
        }
        disp += ((dim)-1)*(int)sqrt(np);
    }
}

/*******************************************************************************
Name:       check_help
Parameters: int argc: number of arguments in argv
            char** argv: command line argument array
            int my_rank: rank of calling process
            int root: root of process communicators
Return:     int: whether the help display is asked for (-1) or not (0)

Desrciption:
Runs through the argument array argv with argc number of elements both given
from commandline to the main function. The function then searches for an
occurrence of the flag '-h'. If it is given the help_Info() function is called
by the root process in order to print the help information on command line and
and an integer -1 is returned. If no help flag is given the function returns 0.
*******************************************************************************/
int check_help(int argc, char **argv, int my_rank, int root){
    int i;
    // catch for help flag
    for(i=0; i<argc; i++){
        if (strcmp( argv[i], "-h") == 0){     // if help flag attached
            /* root process prints information */
            if (my_rank == root){
                help_Info();            // call help information
            }
            return -1;
		}
	}
    return 0;
}

/*******************************************************************************
Name:       help_Info
Parameters:
Return:

Desrciption:
Prints the help information on command line.
*******************************************************************************/
void help_Info(){
    printf("###############################################################\n");
    printf("#                                                             #\n");
    printf("# This is the program to Task 2 of assignment 5               #\n");
    printf("#                                                             #\n");
    printf("# This program can be compiled by following command in its    #\n");
    printf("# encompassing folder:                                        #\n");
    printf("# mpicc -o Task2_master Task2_master.c                        #\n");
    printf("#                                                             #\n");
    printf("# In order to start the program simply write on command line: #\n");
    printf("# one line:                                                   #\n");
    printf("# mpiexec -np 1 ./Task2_master                                #\n");
    printf("#                                                             #\n");
    printf("# Parameters:                                                 #\n");
    printf("#                                                             #\n");
    printf("# -h:       view this help                                    #\n");
    printf("#                                                             #\n");
    printf("###############################################################\n");
}
