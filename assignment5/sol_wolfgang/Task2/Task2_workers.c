/************************************************************************/
/* Program: Task2_workers.c                                             */
/* Author: Wolfgang Stammer <wolfgang.stammer@stud.uni-frankfurt.de>    */
/* matriclenumber: 5962350                                              */
/* Assignment : 5                                                       */
/* Task: 2                                                              */
/* Parameters:  -h: view help                                           */
/*                                                                      */
/* Environment variables: no                                            */
/*                                                                      */
/* Description:                                                         */
/* This program creates an MPI communicator for the worker processes of */
/* Task 2 of Assignment 5. It performs the Cannon Matrix Multiplication */
/* algorithm on the worker processes. Initially every worker process    */
/* receives its submatrix via MPI_Scatterv from the master process, by  */
/* creating a submatrix Datatype. The processes are structured in a     */
/* torus and all communication between processes during the algorithm is*/
/* performed using MPI_Sendrecv_replace and MPI_Cart_shift. Once all    */
/* processes have computed their parts of the result Matrix process 0   */
/* communicates to the master process that the calculations are finished*/
/* and the results can be sent back using MPI_Gatherv.                  */
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

#define TRUE 1
#define FALSE 0
#define WRAP 1

void mult_mat_and_add(double *, double *, double *, int);
void init_ones(int, int[]);
void get_displacements(int, int, int []);

int main(int argc, char* argv[ ]){

    int         my_rank,            // rank of the process
                numprocs,           // number of processes
                dim_SubMat,         // dimension of process submatrix
                dim_Mat,            // dimension of entire matrix
                row_rankA,          // rank in row communicator A
                col_rankB,          // rank in column communicator B
                my_rowID,
                my_colID,
                remain_dims[2],
                rank_src,
                rank_dest;
    double      *my_submatrixA,     // submatrix of matrix A
                *matrixA,           // redundant pointer
                *my_submatrixB,     // submatrix of matrix B
                *matrixB,           // redundant pointer
                *my_submatrixC,     // submatrix of matrix C
                *matrixC;           // redundant pointer
    MPI_Comm    intercomm,
                torus_comm_A,
                torus_row_comm_A,
                torus_comm_B,
                torus_col_comm_B;
    MPI_Status  status;

    /* initializing of MPI-Interface */
    MPI_Init(&argc, &argv);
    /* get rank IDs */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    /* get number of running processes */
    MPI_Comm_size( MPI_COMM_WORLD, &numprocs);
    /* get parent communicator */
    MPI_Comm_get_parent(&intercomm);

    /* get submatrix and matrix dimensions */
    dim_SubMat = strtol(argv[1], NULL, 10); // convert string to int
    dim_Mat = strtol(argv[2], NULL, 10);

    /* allocate memory for all submatrices */
    // matrixA = malloc(dim_Mat * dim_Mat * sizeof(double));
    my_submatrixA = (double*)malloc(dim_SubMat * dim_SubMat * sizeof(double));
    my_submatrixB = (double*)malloc(dim_SubMat * dim_SubMat * sizeof(double));
    my_submatrixC = (double*)calloc(dim_SubMat * dim_SubMat, sizeof(double));

    /* create submatrix datatype */
    MPI_Datatype type, submat_type;
    int sizes[2] = {dim_Mat, dim_Mat},          /* size of global array */
        subsizes[2] = {dim_SubMat, dim_SubMat}, /* size of sub-region */
        starts[2] = {0,0};
    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &type);
    /* change the extent of the type */
    MPI_Type_create_resized(type, 0, dim_SubMat*sizeof(double), &submat_type);
    MPI_Type_commit(&submat_type);

    /* scatter submatrices from root to all worker processes */
    // how many pieces of data everyone has, in units of blocks
    int sndcounts[numprocs];
    init_ones(numprocs, sndcounts);
    // the starting point of everyone's data in the global array, in block extents
    int snddispls[numprocs];
    get_displacements(numprocs, dim_SubMat, snddispls);
    // scatter matrix A
    MPI_Scatterv(matrixA, sndcounts, snddispls, /* proc i gets counts[i] types from displs[i] */
                submat_type, my_submatrixA, dim_SubMat*dim_SubMat, MPI_DOUBLE,
                0, intercomm);
    // scatter matrix B
    MPI_Scatterv(matrixB, sndcounts, snddispls, /* proc i gets counts[i] types from displs[i] */
                submat_type, my_submatrixB, dim_SubMat*dim_SubMat, MPI_DOUBLE,
                0, intercomm);

    /* create matrix A torus */
    int dims = 2,
        dim_size [2] = {sqrt(numprocs), sqrt(numprocs)}, // sqrt(p)xsqrt(p)
        wrap [2] = {TRUE, TRUE},
        reorder = TRUE;
    MPI_Cart_create(MPI_COMM_WORLD, dims, dim_size, wrap, reorder, &torus_comm_A);

    /* create matrix B torus */
    MPI_Cart_create(MPI_COMM_WORLD, dims, dim_size, wrap, reorder, &torus_comm_B);

    /* create row subtopology of matrix A*/
    // first dimension belongs to subtopology, i.e. rows
    remain_dims[0] = FALSE;
    remain_dims[1] = TRUE;
    MPI_Cart_sub(torus_comm_A, remain_dims, &torus_row_comm_A);
    /* get local ranks in row subtopology */
    MPI_Comm_rank(torus_row_comm_A, &row_rankA);

    /* create column subtopology of matrix B */
    // second dimension belongs to subtopology, i.e. columns
    remain_dims[0] = TRUE;
    remain_dims[1] = FALSE;
    MPI_Cart_sub(torus_comm_B, remain_dims, &torus_col_comm_B);
    /* get local ranks in column subtopology */
    MPI_Comm_rank(torus_col_comm_B, &col_rankB);

    /* compute index of my row subtopology of matrix A */
    my_rowID = (int)(my_rank/sqrt(numprocs));
    // printf("rank: %d with row index %d\n", my_rank, my_rowID);

    /* compute index of my column subtopology of matrix B */
    my_colID = my_rank%(int)sqrt(numprocs);
    // printf("rank: %d with column index %d\n", my_rank, my_colID);

    /* calculate source rank and destination rank for initial shift along */
    /* matrix A rows */
    MPI_Cart_shift(torus_row_comm_A, 0, -my_rowID, &rank_src, &rank_dest);
    printf("original rank %d with rank source %d and rank dest %d\n", my_rank, rank_src, rank_dest);

    /* exchange values along rows */
    MPI_Sendrecv_replace(my_submatrixA, dim_SubMat*dim_SubMat, MPI_DOUBLE,
                        rank_dest, 0, rank_src, 0, torus_row_comm_A, &status);

    /* calculate source rank and destination rank for initial shift along */
    /* matrix B columns */
    MPI_Cart_shift(torus_col_comm_B, 0, -my_colID, &rank_src, &rank_dest);
    printf("original rank %d with rank source %d and rank dest %d\n", my_rank, rank_src, rank_dest);

    /* exchange values along columns */
    MPI_Sendrecv_replace(my_submatrixB, dim_SubMat*dim_SubMat, MPI_DOUBLE,
                        rank_dest, 0, rank_src, 0, torus_col_comm_B, &status);

    // make sure all processes have received and sent
    MPI_Barrier(MPI_COMM_WORLD);

    int i, j, k, l, stepIDX;
    /* repeat sqrt(p) steps */
    for(stepIDX = 0; stepIDX < sqrt(numprocs); stepIDX++){
        /* perform actual multiplication */
        mult_mat_and_add(my_submatrixA, my_submatrixB, my_submatrixC, dim_SubMat);
        if(my_rank == 0){
            for (k = 0; k < dim_SubMat; k++){
                for (l = 0; l < dim_SubMat; l++){
                    printf("%.1f, ", my_submatrixC[k * dim_SubMat + l]);
                }
                printf("\n");
            }
        }
        /* shift along row subring */
        MPI_Cart_shift(torus_row_comm_A, 0, -1, &rank_src, &rank_dest);
        MPI_Sendrecv_replace(my_submatrixA, dim_SubMat*dim_SubMat, MPI_DOUBLE,
                            rank_dest, 0, rank_src, 0, torus_row_comm_A, &status);
        /* shift along column subring */
        MPI_Cart_shift(torus_col_comm_B, 0, -1, &rank_src, &rank_dest);
        MPI_Sendrecv_replace(my_submatrixB, dim_SubMat*dim_SubMat, MPI_DOUBLE,
                            rank_dest, 0, rank_src, 0, torus_col_comm_B, &status);
    }

    // make sure all processes have finished calculations
    MPI_Barrier(MPI_COMM_WORLD);

    // if process 0 send to master process that calculations finished
    // and wait for a reply from master process before sending data
    if (my_rank == 0){
//        sleep(5);
        int finished = 1;
        MPI_Request finish_request;

        MPI_Isend(&finished, 1, MPI_INT, 0, 0, intercomm, &finish_request);
        MPI_Wait(&finish_request, &status);
    }

    /* it all goes back to process 0 */
    MPI_Gatherv(
        my_submatrixC,
        dim_SubMat*dim_SubMat,
        MPI_DOUBLE,
        matrixC,
        sndcounts,
        snddispls,
        submat_type,
        0,
        intercomm
    );

    /* free datatype */
    MPI_Type_free(&submat_type);

    /* free pointer array memories */
    free(my_submatrixA);
    free(my_submatrixB);
    free(my_submatrixC);

    MPI_Finalize();		            // finalizing MPI interface

    return 0;
}

/*******************************************************************************
Name:       mult_mat_and_add
Parameters: double *submatrixA: sub matrix of A values
            double *submatrixB: sub matrix of B values
            double *submatrixC: sub matrix of resulting C values
            int dim: dimension of all three square matrices
Return:

Description:
This function computes the matrix multiplication of submatrixA with submatrixB,
i.e. rows times columns and stores the result in submatrixC.
*******************************************************************************/
void mult_mat_and_add(double *submatrixA, double *submatrixB, double *submatrixC, int dim){
    int i, j, k;
    double sum = 0;

    for(i = 0; i < dim; i++){
        for(j = 0; j < dim; j++){
            sum = 0;
            for (k = 0; k < dim; k++){
                sum += submatrixA[i*dim + k] * submatrixB[k*dim + j];
            }
            submatrixC[i*dim + j] += sum;
        }
    }
}

/*******************************************************************************
Name:       init_ones
Parameters: int size: size of array
            int array[size]: array of integers
Return:

Description:
This function receives an integer array and writes it with ones.
*******************************************************************************/
void init_ones(int np, int array[np]){
    int i;
    for (i = 0; i < np; i++){
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
            disp += 1;
        }
        disp += ((dim)-1)*(int)sqrt(np);
    }
}
