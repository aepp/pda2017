/***********************************************************************
Program:            Task2_worker.c
Author:             Joshua Mendaza <joshua_mendaza@outlook.com>
matriclenumber:     4670232
Assignment:         5
Task:               2
Parameters:         no

Environment variables: no

Description:

The program builds an process cartesian topology with other processes.
Each process receives 2 blocks auf matrices A and B and multiplies them.
Then they send block A to the left process and block B to the process above.
after sqrt(process_count) each block is finish.
The result is send to the master process

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
	int t = 1;							// TRUE
    int f = 0;							// FALSE

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



	/// get the intercommunicator of the master process
	// This communicator is used for communication between master and worker processes
	MPI_Comm intercomm;
	MPI_Comm_get_parent(&intercomm);

	/// get metadata from master process
	int metadata[4];				/// receive array
	MPI_Bcast(metadata, 4, MPI_INT, 0, intercomm);

	int matrixLength = metadata[0];	/// total lenth of matrices A and B
    int workerCount = metadata[1];	/// total count of worker
    int workerLength = metadata[2];	/// amount of workers in one row / col
	int secDelay = metadata[3];		/// delay in seconds
    int cellCount = matrixLength * matrixLength;	/// total amount of cells in matrices A and B
    MPI_Status status;				/// communication status


	/// define MPI_Datatype for matrix blocks
	int rows = matrixLength;                /// row count of matrix
    int cols = matrixLength;                /// column count of matrix
    const int procRows = workerLength;      /// number of processes per column
    const int procCols = workerLength;      /// number of processes per row
    const int blockRows = rows / procRows;  /// number of rows in in one block (per process)
    const int blockCols = cols / procCols;  /// number of cols in in one block (per process)
	int blockSize = blockRows * blockCols;  /// number of cells per block

    MPI_Datatype blocktype;					/// datatype of one matrix-block
    MPI_Datatype blocktype2;

    MPI_Type_vector(blockRows, blockCols, cols, MPI_DOUBLE, &blocktype2);
    MPI_Type_create_resized( blocktype2, 0, sizeof(double), &blocktype);
    MPI_Type_commit(&blocktype);

	/// calculate the dicplacements for GatherV
    int disps[workerCount];
    int counts[workerCount];
    for (int ii=0; ii<procRows; ii++) {
        for (int jj=0; jj<procCols; jj++) {
            disps[ii*procCols+jj] = ii*cols*blockRows+jj*blockCols;
            counts [ii*procCols+jj] = 1;
        }
    }


	/// each worker receives one matrix block from the process
	double matrixA[matrixLength][matrixLength];		/// send array, not used by workers
	double matrixB[matrixLength][matrixLength];		/// send array, not used by workers
	double blockA[blockRows][blockCols];			/// block of matrix A, for this process
	double blockB[blockRows][blockCols];			/// block of matrix A, for this process
	MPI_Scatterv(matrixA
                , counts
                , disps
                , blocktype
                , blockA
                , blockRows*blockCols
                , MPI_DOUBLE
                , 0
                , intercomm);

	MPI_Scatterv(matrixB
                , counts
                , disps
                , blocktype
                , blockB
                , blockRows*blockCols
                , MPI_DOUBLE
                , 0
                , intercomm);





	/// create cartesian communicator (grid)
	int my_rank_rows;           // rank of the process in MPI_ROWS
    int my_rank_cols;           // rank of the process in MPI_COLUMNS
    MPI_Comm MPI_CARTESIAN;     // contains the cartesian topology
    MPI_Comm MPI_ROWS;          // contains the row of the cartesion communicator
    MPI_Comm MPI_COLUMNS;       // contains the column of the cartesion communicator
    int dims = 2;               // amount of dimensions
    int dim_size[2] = {workerLength, workerLength};   // size of each dimension (matrix is quadratic)
    int wrap[2] = {1,1};        // no circles allowed
    int reorder = f;            // no reordering
    MPI_Cart_create(MPI_COMM_WORLD, dims, dim_size, wrap, reorder, &MPI_CARTESIAN);

    /// create sub communicator for rows
    int remainDims[2] = {0,1};
    MPI_Cart_sub(MPI_CARTESIAN, remainDims, &MPI_ROWS);
    /// get process-rank for the row communicator
    MPI_Comm_rank(MPI_ROWS, &my_rank_rows);

    /// create sub communicator for columns
    int remainDimsCol[2] = {1,0};
    MPI_Cart_sub(MPI_CARTESIAN, remainDimsCol, &MPI_COLUMNS);
    /// get process-rank for the column communicator
    MPI_Comm_rank(MPI_COLUMNS, &my_rank_cols);



	/// calculate receiver- and sender-rank for initial shift of matrix A
	int rowDisp = -1 * my_rank_cols;
	int rowRankSource;
	int rowRankDest;
	MPI_Cart_shift(MPI_ROWS, 0, rowDisp, &rowRankSource, &rowRankDest);

	/// calculate receiver- and sender-rank for initial shift of matrix B
	int colDisp = -1 * my_rank_rows;
	int colRankSource;
	int colRankDest;
	MPI_Cart_shift(MPI_COLUMNS, 0, colDisp, &colRankSource, &colRankDest);


	/// inital shift of matrix blocks for matrix A
	MPI_Sendrecv_replace(blockA
						, blockSize
						, MPI_DOUBLE
						, rowRankDest
						, 1
						, rowRankSource
						, 1
						, MPI_ROWS
						, &status);

	/// inital shift of matrix blocks for matrix B
	MPI_Sendrecv_replace(blockB
						, blockSize
						, MPI_DOUBLE
						, colRankDest
						, 2
						, colRankSource
						, 2
						, MPI_COLUMNS
						, &status);


	/// declare 2D for block of result matrix C
	double resultBlock[blockRows][blockCols];
	memset(resultBlock, 0.0, sizeof(double) * blockSize);


	/// calculate receiver- and sender-rank for the left/up shift of matrix A and B
	MPI_Cart_shift(MPI_ROWS, 0, -1, &rowRankSource, &rowRankDest);
	MPI_Cart_shift(MPI_COLUMNS, 0, -1, &colRankSource, &colRankDest);


	/// iterate over the square root of the worker count
	for (int iter = 0; iter < workerLength; iter++){

		/// add the matrix-multiplication of block A and block B to block C
		for (int i = 0; i < blockRows; i++){
			for (int j = 0; j < blockRows; j++){
				for (int k = 0; k < blockRows; k++){
					resultBlock[i][j] += blockA[i][k] * blockB[k][j];
				}
			}
		}

		/// Send and receive the block of matrix A
		MPI_Sendrecv_replace(blockA
							, blockSize
							, MPI_DOUBLE
							, rowRankDest
							, 1
							, rowRankSource
							, 1
							, MPI_ROWS
							, &status);


		/// Send and receive the block of matrix B
		MPI_Sendrecv_replace(blockB
							, blockSize
							, MPI_DOUBLE
							, colRankDest
							, 2
							, colRankSource
							, 2
							, MPI_COLUMNS
							, &status);

	}

	/// Wait for all processes to finish
	MPI_Barrier(MPI_COMM_WORLD);

	/// Process 0 is used to send finish-status to the master process
	if (my_rank == 0){
		/// delay the worker-group (for testing)
		sleep(secDelay);

		/// start a non-blocking communication with the master process
		/// the master process can use MPI_Test for this receive
		///, to check if this worker-group is already finished
		int finished = 1;
		MPI_Request sendFinished;                 // request handles for pers. commmunication
		MPI_Isend(&finished, 1, MPI_INT, 0, 1, intercomm, &sendFinished);
		MPI_Wait(&sendFinished, &status);
	}

	/// Now each process sends its result-block of matrix C to the master process
	double resultC[matrixLength][matrixLength];
    MPI_Gatherv(resultBlock, blockSize, MPI_DOUBLE, resultC, counts, disps, blocktype, 0, intercomm);

    MPI_Finalize();
    return 0;
}
