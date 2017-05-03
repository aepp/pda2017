/************************************************************************/
/* Program: assignment.c                                                */
/* Author: Aleksandr Epp <aleksandr.epp@gmail.com>                      */
/* Matriclenumber: 6002853                                              */
/* Assignment : 1                                                       */
/* Parameters: -h, -t, -m                                               */
/* Environment variables: no                                            */
/*                                                                      */
/* Description:                                                         */
/*                                                                      */
/* This is the main program of Assignment 1. From here individual tasks */
/* can be started. Or view the help.                                    */
/*                                                                      */
/************************************************************************/ 

#include <stdio.h>      // import of the definitions of the C IO library
#include <stdlib.h>     // to use exit() function
#include <string.h>     // import of the definitions of the string operations
#include <unistd.h>     // standard unix io library definitions and declarations
#include <errno.h>      // system error numbers
#include <inttypes.h>   // for strtoumax() string to int coversion

#include "help.h"       // to show help with -h parameter
#include "task1.h"      // to run Task 1
#include "task2a.h"     // to run Task 2a)
#include "task2b.h"     // to run Task 2b)

// constants
#define DEFAULT_ARRAY_SIZE 100 // default number of elements in random int array
#define MAX_RANDOM 100 // maximum random integer that appears in the array
#define ASSIGNMENT_NR 1

int main(int argc, char* argv[])
{ 
    int arraySize = DEFAULT_ARRAY_SIZE, // number of elements in random int array
        opt = 0; // command line parameter name

    char *taskName = NULL; // task nr. given by cli parameter

    // iterate through all cli parameters
    while ((opt = getopt(argc, argv, "m:t:h")) != -1) {
        switch(opt) {
            case 'h': // if -h parameter given => show help
                printf("Show help:\n");
                showHelp(ASSIGNMENT_NR);
                exit(0);
            case 'm': // if -m parameter given => set arraySize
                arraySize = strtoumax(optarg, NULL, 10);
                break;
            case 't': // if -t parameter given => save task nr. to decide which task to execute
                taskName = optarg;
                break;
            default: // if any other cli parameters provided show
                // exit program with error status
                exit(1);
        }
    }

    // decide which task to run based on -t parameter
    if(taskName == NULL){ // if no parameter provided - just show help
        showHelp(ASSIGNMENT_NR);
        exit(0);
    } else if(!strcmp(taskName, "1")) {
        task1(&argc, &argv);
    } else if(!strcmp(taskName, "2a")) {
        task2a(argc, argv, MAX_RANDOM, arraySize);
    } else if(!strcmp(taskName, "2b")) {
        task2b(argc, argv, MAX_RANDOM, arraySize);
    } else { // if invalid task number provided
        printf("Task %s doesn't exist. See help (-h) for available tasks.\n", taskName);
    }

    // end of program with exit code 0
    return 0;
}