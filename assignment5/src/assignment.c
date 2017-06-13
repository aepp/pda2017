/************************************************************************/
/* Program: assignment.c                                                */
/* Author: Aleksandr Epp <aleksandr.epp@gmail.com>                      */
/* Matriclenumber: 6002853                                              */
/* Assignment : 5                                                       */
/* Parameters: -m, -h                       see help (-h) for details   */
/* Environment variables: no                                            */
/*                                                                      */
/* Description:                                                         */
/*                                                                      */
/* This is the main program of Assignment 5. From here individual tasks */
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

// constants & default values
#define ASSIGNMENT_NR       5
#define RANDOM_INT_COUNT    1           // default amount of random numbers per process

int main(int argc, char* argv[])
{
    int opt = 0,                        // command line parameter name
        randomIntCount = RANDOM_INT_COUNT;    // amount of random numbers per process

    // iterate through all cli parameters
    while ((opt = getopt(argc, argv, "m:h")) != -1) {
        switch(opt) {
            case 'h': // if -h parameter given => show help
                showHelp(ASSIGNMENT_NR, argc, argv);
                exit(0);
            case 'm': // if -m parameter given => set the amount of random numbers per process
                randomIntCount = strtoumax(optarg, NULL, 10);
                break;
            default: // if any other cli parameters provided
                // exit program with error status
                exit(1);
        }
    }

    task1(argc, argv, randomIntCount);

    // end of program with exit code 0
    return 0;
}