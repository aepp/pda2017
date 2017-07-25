/************************************************************************/
/* Program: Task1.c                                                     */
/* Author: Wolfgang Stammer <wolfgang.stammer@stud.uni-frankfurt.de>    */
/* matriclenumber: 5962350                                              */
/* Assignment : 5                                                       */
/* Task: 1                                                              */
/* Parameters:  -h: view help                                           */
/*                                                                      */
/* Environment variables: no                                            */
/*                                                                      */
/* Description:                                                         */
/* This program applies a specific filter to a 2D input image and saves */
/* the resulting image to a specified output file. Every process reads  */
/* an equal amount of rows from the original image and converts the     */
/* initially unsigned char values to double. The individual arrays are  */
/* then padded by two columns and two rows on each side of the 2D       */
/* arrays. The newly padded first two rows are the last two rows of the */
/* original image array that the previous neighboring process has,      */
/* whereas the newly padded last two rows are the first two rows of the */
/* original image array that the following neighboring process has.     */
/* This sending of rows is done using persistent communication.         */
/* These padded arrays are then convolved with a 5x5 filter array, the  */
/* type of which is specified by the user (e.g. blur, relief, ...).     */
/* The final filtered image is converted back to unsigned char and      */
/* in the output file.                                                  */
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

#define sizeFilter 5
#define padSize 4

double apply_Filter_once(double *, double *, int, int, int);
void apply_Filter(double *, double *, double *, int, int);
void pad_Array_top_bottom(double *, double *, double *, double *, int, int);
void get_Filter(double *, char[]);
void conv_CharArray2DoubleArray(unsigned char*, double*, int, int);
void conv_DoubleArray2UnsignCharArray(unsigned char*, double*, int, int);
void get_strParam(int, char**, char[], char[]);
int get_intParam(int, char**, char[]);
int check_help(int, char**, int, int);
void help_Info();

int main(int argc, char* argv[ ]){

    const int root = 0;
    int     my_rank, 					// rank of the process
            numprocs,                   // number of processes
            err = 0,
            dim1,                       // number of rows in image
            dim2,                       // number of columns in image
            dim1PerProc,                // number of rows per process
            numIters,                // number of filter iterations (strength)
            i, j;
    double  tmp_dim1PerProc,            // number of rows per process
            *my_RowsDouble,             // image row values in double
            *my_RowsDoublePadded,       // result values in double
            *filterArray,               // filter array
            *snd_topArray,              // array for top rows pers. comm.
            *snd_bottomArray,           // array for bottom rows pers. comm.
            *rcv_topArray,              // array for top rows pers. comm.
            *rcv_bottomArray,           // array for bottom rows pers. comm.
            *dummy_array;               // dummy array of zeros
    char    inputFileName[PATH_MAX],            // file path for Image
            outputFileName[PATH_MAX],           // output file Name
            filter_type[PATH_MAX];
    unsigned char   *my_RowsChar;       // original image values in unsign char
    MPI_File fh;                        // file handle
    MPI_Status status;
    MPI_Offset fileSize;

    /* set memory for file path variables */
    memset(inputFileName, 0, PATH_MAX);
    memset(outputFileName, 0, PATH_MAX);

    /* initializing of MPI-Interface */
	MPI_Init(&argc, &argv);
    /* get rank IDs */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    /* get number of running processes */
    MPI_Comm_size( MPI_COMM_WORLD, &numprocs);

    /* check for help request and if necessary finalize */
    if (check_help(argc, argv, my_rank, root) == -1){
        MPI_Finalize();
        return -1;
    }

    /* get number of pixels per row */
    // dim2 = 1280;
    dim2 = get_intParam(argc, argv, "-dim2");

    /* get number of iterations */
    numIters = get_intParam(argc, argv, "-m");
    // numIters = 2;

    /* Get file path for image file and output file */
    // strcpy(inputFileName, "../5/ffm_1280x960.gray");
    get_strParam(argc, argv, "-fin", inputFileName);
    // strcpy(outputFileName, "Out.gray");
    get_strParam(argc, argv, "-fout", outputFileName);

    /* get filter type */
    // strcpy(filter_type, "relief");
    get_strParam(argc, argv, "-ftype", filter_type);

    if (my_rank == root) {
        printf("\nThe file path to the input image: %s\n", inputFileName);
        printf("The file path to the output image: %s\n", outputFileName);
        printf("The filter strength: %d \n", numIters);
        printf("Finally the specified filter: %s \n \n", filter_type);
    }

    /* open and read image from file using ordered reading */
    err = MPI_File_open(MPI_COMM_WORLD, inputFileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if (err) MPI_Abort(MPI_COMM_WORLD, 1);      // exit on error

    /* get size of file */
    MPI_File_get_size(fh, &fileSize);
    fileSize = fileSize/sizeof(unsigned char);  // number of pixels altogether
    dim1 = fileSize/dim2;          // number or rows altogether
    tmp_dim1PerProc = (double)dim1/numprocs;  // number or rows per process

    /* check if number of rows can be divided by number of processes */
    if ((int)tmp_dim1PerProc != tmp_dim1PerProc){
        // root process prints message if number of processes incorrect
        if (my_rank == root){
            printf("\nThe number of rows must be equally dividable by the number of processes. \n"
                    "The number of rows for this image is %d. Please restart the program \n"
                    "with a different number of processes. \n \n", dim1);
        }
        MPI_File_close(&fh);                        // close file
        MPI_Finalize();
        return -2;
    }
    // else cast to integer
    else dim1PerProc = (int)tmp_dim1PerProc;

    /* allocate filterArray and get specific predefined filterarray */
    filterArray = malloc(sizeFilter * sizeFilter *sizeof(double));
    get_Filter(filterArray, filter_type);

    /* read rows from file */
    my_RowsChar = malloc(dim1PerProc * dim2 * sizeof(unsigned char));
    // collective blocking
    err = MPI_File_read_ordered(fh, my_RowsChar, dim1PerProc*dim2, MPI_UNSIGNED_CHAR, &status);
    if (err) MPI_Abort(MPI_COMM_WORLD, 1);      // exit on error
    MPI_File_close(&fh);                        // close file

    /* allocate double array and zero padded array */
    my_RowsDouble = malloc(dim1PerProc * dim2 * sizeof(double));
    my_RowsDoublePadded = malloc((dim1PerProc+padSize) * (dim2+padSize) * sizeof(double));
    /* convert char to double */
    conv_CharArray2DoubleArray(my_RowsChar, my_RowsDouble, dim1PerProc, dim2);

    /* initialize requests for persistent communication */
    MPI_Request snd_top_req, snd_bottom_req, rcv_top_req, rcv_bottom_req;

    /* allocate memory for send and receive arrays */
    snd_topArray = malloc(2*dim2*sizeof(double));
    snd_bottomArray = malloc(2*dim2*sizeof(double));
    rcv_topArray = malloc(2*dim2*sizeof(double));
    rcv_bottomArray = malloc(2*dim2*sizeof(double));
    /* allocate dummy array of zeros */
    dummy_array = calloc(2*dim2, sizeof(double));

    // in case of only one process there is no communication
    if (numprocs > 1){
        /* initialize persistent communication */
        /* the order of sending and receiving: the top rows of my neighbour is my */
        /* bottom row, the bottom row of my neighbour is my top row */
        if (my_rank == root){
            // send top dummy rows (zeros) to last process
            MPI_Send_init(dummy_array, 2 * dim2, MPI_DOUBLE,  numprocs-1, 0, MPI_COMM_WORLD, &snd_top_req);
            MPI_Recv_init(rcv_topArray, 2 * dim2, MPI_DOUBLE, numprocs-1, 0, MPI_COMM_WORLD, &rcv_top_req);
            MPI_Send_init(snd_bottomArray, 2 * dim2, MPI_DOUBLE,  my_rank+1, 0, MPI_COMM_WORLD, &snd_bottom_req);
            MPI_Recv_init(rcv_bottomArray, 2 * dim2, MPI_DOUBLE, my_rank+1, 0, MPI_COMM_WORLD, &rcv_bottom_req);
        }
        else if (my_rank == (numprocs-1)){      // last process
            // send bottom dummy rows (zeros) to first process
            MPI_Send_init(snd_topArray, 2 * dim2, MPI_DOUBLE,  my_rank-1, 0, MPI_COMM_WORLD, &snd_top_req);
            MPI_Recv_init(rcv_topArray, 2 * dim2, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD, &rcv_top_req);
            MPI_Send_init(dummy_array, 2 * dim2, MPI_DOUBLE,  root, 0, MPI_COMM_WORLD, &snd_bottom_req);
            MPI_Recv_init(rcv_bottomArray, 2 * dim2, MPI_DOUBLE, root, 0, MPI_COMM_WORLD, &rcv_bottom_req);
        }
        else {
            MPI_Send_init(snd_topArray, 2 * dim2, MPI_DOUBLE,  my_rank-1, 0, MPI_COMM_WORLD, &snd_top_req);
            MPI_Recv_init(rcv_topArray, 2 * dim2, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD, &rcv_top_req);
            MPI_Send_init(snd_bottomArray, 2 * dim2, MPI_DOUBLE,  my_rank+1, 0, MPI_COMM_WORLD, &snd_bottom_req);
            MPI_Recv_init(rcv_bottomArray, 2 * dim2, MPI_DOUBLE, my_rank+1, 0, MPI_COMM_WORLD, &rcv_bottom_req);
        }
    }

    /* apply filter numIters times for the strength of the filter */
    for (k = 0; k < numIters; k++){

        // start and wait for communication
        if (numprocs > 1){
            /* copy rows from my_RowsDouble to send communication arrays */
            memcpy(snd_topArray, my_RowsDouble, 2*dim2*sizeof(double));
            memcpy(snd_bottomArray, my_RowsDouble+(dim1PerProc-2)*dim2, 2*dim2*sizeof(double));

            /* send and receive overlapping image rows with persistent communication */
            // start communication
            MPI_Start(&snd_top_req);
            MPI_Start(&snd_bottom_req);
            MPI_Start(&rcv_top_req);
            MPI_Start(&rcv_bottom_req);
            // wait for successful communication
            MPI_Wait(&snd_top_req, &status);
            MPI_Wait(&snd_bottom_req, &status);
            MPI_Wait(&rcv_top_req, &status);
            MPI_Wait(&rcv_bottom_req, &status);
        }
        // in case of only one process there is no communication
        else{
            memcpy(rcv_topArray, dummy_array, 2*dim2*sizeof(double));
            memcpy(rcv_bottomArray, dummy_array, 2*dim2*sizeof(double));
        }

        /* pad array with new specified top and bottom arrays */
        pad_Array_top_bottom(my_RowsDouble, my_RowsDoublePadded, rcv_topArray, rcv_bottomArray, \
                            dim1PerProc+padSize, dim2+padSize);

        /* apply filter convolution over entire array once */
        apply_Filter(my_RowsDouble, my_RowsDoublePadded, filterArray, dim1PerProc, dim2);
    }

    /* convert double result array back to unsigned char */
    conv_DoubleArray2UnsignCharArray(my_RowsChar, my_RowsDouble, dim1PerProc, dim2);


    /* write result arrays to output file */
    err = MPI_File_open(MPI_COMM_WORLD, outputFileName, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    if (err) MPI_Abort(MPI_COMM_WORLD, 1);      // exit on error
    err = MPI_File_write_ordered(fh, my_RowsChar, dim1PerProc*dim2, MPI_UNSIGNED_CHAR, &status);
    if (err) MPI_Abort(MPI_COMM_WORLD, 1);      // exit on error
    MPI_File_close(&fh);                        // close file


    if (my_rank == root){
        printf("The filtered image has been saved in the specified file %s.\n \n"
                "In order to see the image use the command: \n"
                "display -depth 8 -size 1280xNNNN %s \n \n", outputFileName, outputFileName);
    }

    /* free requests of persistent communication */
    if (numprocs > 1){
        MPI_Request_free(&snd_top_req);
        MPI_Request_free(&snd_bottom_req);
        MPI_Request_free(&rcv_top_req);
        MPI_Request_free(&rcv_bottom_req);
    }

    /* free memory of pointer arrays */
    free(snd_topArray);
    free(snd_bottomArray);
    free(rcv_topArray);
    free(rcv_bottomArray);
    free(filterArray);
    free(my_RowsChar);
    free(my_RowsDouble);
    free(my_RowsDoublePadded);

    MPI_Finalize();		            // finalizing MPI interface

    return 0;
}

/*******************************************************************************
Name:       apply_Filter_once
Parameters: double *my_RowsDoublePadded
            double *filterArray
            int IDX1
            int IDX2
            int dim2
Return:

Description:
This function receives a array, my_RowsDoublePadded, which is padded by two
columns of zeros in the first two and last two columns and specific, possibly
non-zero values in the first two and last two rows. The function then computes
the convolution of the filterArray with my_RowsDoublePadded at the index
[IDX1, IDX2]. The filterArray is then layed over my_RowsDoublePadded, where
[IDX1, IDX2] is at the center of the filter. dim2 specifies the number of
columns in my_RowsDoublePadded.
*******************************************************************************/
double apply_Filter_once(double *my_RowsDoublePadded, double *filterArray,
                            int IDX1, int IDX2, int dim2){
    int v, u,
        shiftIDX1,
        shiftIDX2;
    double result = 0;

    for (v = 0; v < sizeFilter; v++){
        for (u = 0; u < sizeFilter; u++){
            /* shifted indices of padded array */
            shiftIDX1 = IDX1 + v - 2;
            shiftIDX2 = IDX2 + u - 2;

            result += filterArray[v * sizeFilter + u] * my_RowsDoublePadded[shiftIDX1 * dim2 + shiftIDX2];
        }
    }
    return result;
}

/*******************************************************************************
Name:       apply_Filter
Parameters: double *my_RowsDouble
            double *my_RowsDoublePadded
            double *filterArray
            int dim1
            int dim2
Return:

Description:
This function iterates over every element of the array my_RowsDouble and
convolves the corresponding value in the padded array, my_RowsDoublePadded, with
the filter array, filterArray. In this way every value of my_RowsDoublePadded
that is in the non padded area is at the center of the filter array once. The
actual convolution is performed with the function apply_Filter_once().
If the resulting value is too large or too small it is limited and the final
value is stored back into my_RowsDouble.
*******************************************************************************/
void apply_Filter(double *my_RowsDouble, double *my_RowsDoublePadded,
                    double *filterArray, int dim1, int dim2){
    int i, j;
    double tmp;

    /* go through non padded area of paddedArray and apply filter with */
    /* current pixel as center of filter */
    for (i = 0; i < dim1; i++){
        for (j = 0; j < dim2; j++){
            // apply the filter across all pixels of the process once
            // the indices are shifted by two to correct for the different
            // dimensions of the padded array
            // the padding size is added to padded array dimension
            tmp = apply_Filter_once(my_RowsDoublePadded, filterArray, i+2, j+2, dim2+padSize);

            // if outside of gray boundaries (i.e. 0 > tmp or tmp > 255)
            // set accordingly
            if (tmp < 0) tmp = 0;
            else if (tmp > 255) tmp = 255;

            my_RowsDouble[i*dim2 + j] = tmp;
        }
    }
}

/*******************************************************************************
Name:       pad_Array_top_bottom
Parameters: double *origArray
            double *paddedArray
            double *topArray
            double *bottomArray
            int dim1
            int dim2
Return:

Description:
This function receives a 2D array, origArray, and writes its values into a new
array, paddedArray, in the subarray that is not the first and last two rows and
columns. The first and last two columns of paddedArray are set to zeros. Lastly
the columns 2 to (dim2-2) of the the first two rows are set to the values of
topArray and the same columns of the last two rows are set to the values of
bottomArray. dim1 and dim2 specify the number of rows and columns, respectively.
*******************************************************************************/
void pad_Array_top_bottom(double *origArray, double *paddedArray, double *topArray, \
                        double *bottomArray, int dim1, int dim2){
    int i,j;

    for(i = 0; i < dim1; i++){
        for(j = 0; j < dim2; j++){

            // if in first two rows copy from topArray
            if (i < 2 & j >= 2 & j < (dim2-2)){
                paddedArray[i*dim2 + j] = topArray[i*(dim2-padSize) + (j-2)];
            }
            // if in last two rows copy from bottomArray
            else if (i >= (dim1-2) & j >= 2 & j < (dim2-2)){
                paddedArray[i*dim2 + j] = bottomArray[(i-(dim1-2))*(dim2-padSize) + (j-2)];
            }
            // if indices are within the range of the original region
            else if(i >= 2 & i <(dim1-2) & j >= 2 & j <(dim2-2)){
                paddedArray[i*dim2 + j] = origArray[(i-2)*(dim2-padSize) + (j-2)];
            }
            // otherwise pad with zero
            else {
                paddedArray[i*dim2 + j] = 0;
            }
        }
    }
}


/*******************************************************************************
Name:       get_Filter
Parameters: double *filterArray
            char filter_type[]
Return:

Description:
This function receives an array pointer, filterArray, and a string, filter_type.
Several different 2D filter arrays are defined. Depending on which filter is
requested, specified with filter_type, the corresponding filter array is copied
to the array pointer filterArray.
*******************************************************************************/
void get_Filter(double *filterArray, char filter_type[]){

    /* define several different possible 2D filter arrays */
    static const double blurArray[sizeFilter][sizeFilter] = {
                                            {0, 0, 1./37, 0, 0},
                                            {0, 2./37, 4./37, 2./37, 0},
                                            {1./37, 4./37, 9./37, 4./37, 1./37},
                                            {0, 2./37, 4./37, 2./37, 0},
                                            {0, 0, 1./37, 0, 0}};
    static const double sharpenArray[sizeFilter][sizeFilter] = {
                                            {0, 0, 0, 0, 0},
                                            {0, 0, -1, 0, 0},
                                            {0, -1, 5, -1, 0},
                                            {0, 0, -1, 0, 0},
                                            {0, 0, 0, 0, 0}};
    static const double reliefArray[sizeFilter][sizeFilter] = {
                                            {0, 0, 0, 0, 0},
                                            {0, -2, -1, 0, 0},
                                            {0, -1, 1, 1, 0},
                                            {0, 0, 1, 2, 0},
                                            {0, 0, 0, 0, 0}};
    static const double edgeArray[sizeFilter][sizeFilter] = {
                                            {0, 0, 0, 0, 0},
                                            {0, 1./4, 2./4, 1./4, 0},
                                            {0, 2./4, -12./4, 2./4, 0},
                                            {0, 1./4, 2./4, 1./4, 0},
                                            {0, 0, 0, 0, 0}};
    static const double testArray[sizeFilter][sizeFilter] = {
                                            {0, 0, 0, 0, 0},
                                            {0, 0, 0, 0, 0},
                                            {0, 0, 1, 0, 0},
                                            {0, 0, 0, 0, 0},
                                            {0, 0, 0, 0, 0}};


    /* depending on which filter is asked for, copy local array (filter) to
    filterArray pointer */
    if(strcmp(filter_type,"blur") == 0){
        memcpy(filterArray, blurArray, sizeof blurArray);
    }
    else if(strcmp(filter_type,"sharpen") == 0){
        memcpy(filterArray, sharpenArray, sizeof sharpenArray);
    }
    else if(strcmp(filter_type,"relief") == 0){
        memcpy(filterArray, reliefArray, sizeof reliefArray);
    }
    else if(strcmp(filter_type,"edge") == 0){
        memcpy(filterArray, edgeArray, sizeof edgeArray);
    }
    else if(strcmp(filter_type,"test") == 0){
        memcpy(filterArray, testArray, sizeof testArray);
    }
}

/*******************************************************************************
Name:       conv_CharArray2DoubleArray
Parameters: unsigned char *charArray
            double *doubleArray
            int dim1
            int dim2
Return:

Description:
This function receives an unsigned chararray and iterates through all values,
casting them as doubles and saving these in the array doubleArray.
*******************************************************************************/
void conv_CharArray2DoubleArray(unsigned char *charArray, double *doubleArray,
                                int dim1, int dim2){
    int i,j;

    for (i = 0; i < dim1; i++){
        for (j = 0; j < dim2; j++){
            doubleArray[i * dim2 + j] = (double)charArray[i * dim2 + j];
        }
    }
}

/*******************************************************************************
Name:       conv_DoubleArray2UnsignCharArray
Parameters: unsigned char *charArray
            double *doubleArray
            int dim1
            int dim2
Return:

Description:
This function receives a double array and iterates through all values, casting
them as unsigned chars and saving these in the array charArray.
*******************************************************************************/
void conv_DoubleArray2UnsignCharArray(unsigned char *charArray, double *doubleArray,
                                int dim1, int dim2){
    int i,j;

    for (i = 0; i < dim1; i++){
        for (j = 0; j < dim2; j++){
            charArray[i * dim2 + j] = (unsigned char)doubleArray[i * dim2 + j];
        }
    }
}

/*******************************************************************************
Name:       get_strParam
Parameters: int argc
            char** argv
            char param[]
            char returnStr[]
Return:

Description:
Runs through the argument array argv with argc number of elements both given
from commandline to the main function. The function then searches for an
occurrence of the flag given in param[]. The next parameter is the relevant
value given on commandline and is copied to returnStr[].
*******************************************************************************/
void get_strParam(int argc, char **argv, char param[], char returnStr[]) {
    int i;
    char c;

    for(int i=0; i<argc; i++){  // go through argument variables to catch m flag
        if (strcmp( argv[i], param) == 0){
            strcpy(returnStr, argv[i+1]);
            break;
        }
    }
}

/*******************************************************************************
Name:       get_intParam
Parameters: int argc
            char** argv
            char param[]
Return:     int returnInt

Description:
Runs through the argument array argv with argc number of elements both given
from commandline to the main function. The function then searches for an
occurrence of the flag given in param[]. The next parameter is the relevant
value given on commandline and is returned.
*******************************************************************************/
int get_intParam(int argc, char **argv, char param[]) {
    int i, returnInt;

    for(int i=0; i<argc; i++){  // go through argument variables to catch m flag
        if (strcmp( argv[i], param) == 0){
            sscanf (argv[i+1],"%d",&returnInt);         // convert character to int
            break;
        }
    }
    return returnInt;
}

/*******************************************************************************
Name:       check_help
Parameters: int argc
            char** argv
            int my_rank
            int root
Return:     int

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
    printf("# This is the program to Task 1 of assignment 5               #\n");
    printf("#                                                             #\n");
    printf("# This program can be compiled by following command in its    #\n");
    printf("# encompassing folder:                                        #\n");
    printf("# mpicc -o Task1 Task1.c                                      #\n");
    printf("#                                                             #\n");
    printf("# In order to start on all available CPUs, simply write in    #\n");
    printf("# one line:                                                   #\n");
    printf("# mpiexec ./Task1 -fin fin -fout fout -m m -dim2 dim2         #\n");
    printf("# -ftype ftype                                                #\n");
    printf("# The parameters fin, fout, m, dim2 and ftype are described   #\n");
    printf("# below.                                                      #\n");
    printf("#                                                             #\n");
    printf("# In order to start the program for a specific number of CPUs #\n");
    printf("# use following command, where np specifies the number:       #\n");
    printf("# mpiexec -np np ./Task1 -fin fin -fout fout -m m -dim2 dim2  #\n");
    printf("# -ftype ftype                                                #\n");
    printf("# Please use a square number of processes, i.e. np = nxn,     #\n");
    printf("# e.g. 9.                                                     #\n");
    printf("#                                                             #\n");
    printf("# Parameters:                                                 #\n");
    printf("#                                                             #\n");
    printf("# -h:       view this help                                    #\n");
    printf("# -fin:     path to input file                                #\n");
    printf("# -fout:    path to output file                               #\n");
    printf("# -m:       integer, number of filter iterations              #\n");
    printf("#           (strength of filter)                              #\n");
    printf("# -dim2:    integer, number of elements along 2nd dimension   #\n");
    printf("#           (columns) of input image                          #\n");
    printf("# -ftype:   filter type to be used, please use one of the     #\n");
    printf("#           following:                                        #\n");
    printf("#           - blur                                            #\n");
    printf("#           - sharpen                                         #\n");
    printf("#           - relief                                          #\n");
    printf("#           - edge                                            #\n");
    printf("#                                                             #\n");
    printf("###############################################################\n");
}
