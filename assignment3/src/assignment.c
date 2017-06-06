/************************************************************************/
/* Program: assignment.c                                                */
/* Author: Aleksandr Epp <aleksandr.epp@gmail.com>                      */
/* Matriclenumber: 6002853                                              */
/* Assignment : 3                                                       */
/* Parameters: -h, -t, -m, -v, -s, -e       see help (-h) for details   */
/* Environment variables: no                                            */
/*                                                                      */
/* Description:                                                         */
/*                                                                      */
/* This is the main program of Assignment 3. From here individual tasks */
/* can be started. Or view the help.                                    */
/*                                                                      */
/************************************************************************/

#include <stdio.h>                      // import of the definitions of the C IO library
#include <stdlib.h>                     // to use exit() function
#include <string.h>                     // import of the definitions of the string operations
#include <unistd.h>                     // standard unix io library definitions and declarations
#include <errno.h>                      // system error numbers
#include <inttypes.h>                   // for strtoumax() string to int conversion

#include "help.h"                       // to show help with -h parameter
#include "task1.h"                      // to run Task 1
#include "task2.h"                      // to run Task 2

#include <linux/limits.h>               // for PATH_MAX

// constants & default values
#define ASSIGNMENT_NR   3
#define EPSILON         0.000001       // default epsilon value for Jacobi method
#define MATRIX          ((unsigned char *)"./examples/Matrix_A_8x8")    // default matrix A path
#define VECTOR_B        ((unsigned char *)"./examples/Vector_b_8x")     // default vector b path
#define VECTOR_B_SIZE   8               // default vector b size

int main(int argc, char* argv[])
{
    int opt = 0,                        // command line parameter name
        vectorBSize = VECTOR_B_SIZE;    // size of vector b
    double epsilon = EPSILON;           // epsilon value for Jacobi method

    char *taskName = NULL,              // task nr. given by cli parameter
         *matrixAFilePath = MATRIX,     // path to the file containing matrix A
         *vectorBFilePath = VECTOR_B;   // path to the file containing vector b

//    memset(matrixAFilePath, 0, PATH_MAX);
//    memset(vectorBFilePath, 0, PATH_MAX);
//
//    strcpy(matrixAFilePath, "./examples/Matrix_A_8x8");
//    strcpy(vectorBFilePath, "./examples/Vector_b_8x");

    // iterate through all cli parameters
    while ((opt = getopt(argc, argv, "t:e:m:v:s:h")) != -1) {
        switch(opt) {
            case 'h': // if -h parameter given => show help
                showHelp(ASSIGNMENT_NR, EPSILON, argc, argv);
                exit(0);
            case 't': // if -t parameter given => save task nr. to decide which task to execute
                taskName = optarg;
                break;
            case 'm': // if -m parameter given => set path to matrix file
                matrixAFilePath = optarg;
//                strcpy(matrixAFilePath, optarg);
                break;
            case 'v': // if -v parameter given => set path to vector b file
                vectorBFilePath = optarg;
//                strcpy(vectorBFilePath, optarg);
                break;
            case 's': // if -s parameter given => set the size of vector b
                vectorBSize = strtoumax(optarg, NULL, 10);
                break;
            case 'e': // if -e parameter given => set epsilon
                sscanf(optarg, "%lf", &epsilon);
                break;
            default: // if any other cli parameters provided
                // exit program with error status
                exit(1);
        }
    }

    // decide which task to run based on -t parameter
    if(taskName == NULL){ // if no parameter provided - just show help
        showHelp(ASSIGNMENT_NR, EPSILON, argc, argv);
        exit(0);
    } else if(!strcmp(taskName, "1")) {
        task1(argc, argv, epsilon, matrixAFilePath, vectorBFilePath, vectorBSize);
    } else if(!strcmp(taskName, "2")) {
        task2(argc, argv, epsilon, matrixAFilePath, vectorBFilePath, vectorBSize);
    } else { // if invalid task number provided
        printf("Task %s doesn't exist. See help (-h) for available tasks.\n", taskName);
    }

    // end of program with exit code 0
    return 0;
}