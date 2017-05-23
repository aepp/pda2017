/************************************************************************/
/* Program: assignment.c                                                */
/* Author: Aleksandr Epp <aleksandr.epp@gmail.com>                      */
/* Matriclenumber: 6002853                                              */
/* Assignment : 1                                                       */
/* Parameters: -h, -t, -m, -n, -a ,-b, -f                               */
/* Environment variables: no                                            */
/*                                                                      */
/* Description:                                                         */
/*                                                                      */
/* This is the main program of Assignment 2. From here individual tasks */
/* can be started. Or view the help.                                    */
/*                                                                      */
/************************************************************************/ 

#include <stdio.h>          // import of the definitions of the C IO library
#include <stdlib.h>         // to use exit() function
#include <string.h>         // import of the definitions of the string operations
#include <unistd.h>         // standard unix io library definitions and declarations
#include <errno.h>          // system error numbers
#include <inttypes.h>       // for strtoumax() string to int conversion

#include "help.h"           // to show help with -h parameter
#include "task1.h"          // to run Task 1
#include "task2.h"          // to run Task 2

// constants
#define DEFAULT_SIZE 10    // default size of random int array or amount of intervals
#define ASSIGNMENT_NR 2

int main(int argc, char* argv[])
{ 
    int opt = 0; // command line parameter name

    double a = 0, b = 0, // integration limits
           n = DEFAULT_SIZE, // default n for both tasks
           func = 1, // default function to integrate for task 1
           commMode = 1; // default communication mode for task 2

    char *taskName = NULL; // task nr. given by cli parameter

    // iterate through all cli parameters
    while ((opt = getopt(argc, argv, "hn:t:a:b:f:m:")) != -1) {
        switch(opt) {
            case 'h': // if -h parameter given => show help
                printf("Show help:\n");
                showHelp(ASSIGNMENT_NR);
                exit(0);
            case 'n': 
                // if -n parameter given => 
                //      set number of intervals per process for task 1
                //      or set number of random integers for task 2
                sscanf(optarg, "%lf", &n);
                break;
            case 't': // if -t parameter given => save task nr. to decide which task to execute
                taskName = optarg;
                break;
            case 'a': // if -a parameter given => set left integration limit
                sscanf(optarg, "%lf", &a);
                break;
            case 'b': // if -b parameter given => set right integration limit
                sscanf(optarg, "%lf", &b);
                break;
            case 'f': // if -f parameter given => select which function to integrate
                sscanf(optarg, "%lf", &func);
                break;
            case 'm': // if -f parameter given => select communication mode for task 2
                sscanf(optarg, "%lf", &func);
                break;
            default: // if any other cli parameters provided
                // exit program with error status
                exit(1);
        }
    }

    // decide which task to run based on -t parameter
    if(taskName == NULL){ // if no parameter provided - just show help
        showHelp(ASSIGNMENT_NR);
        exit(0);
    } else if(!strcmp(taskName, "1")) {
        if(b != 0){
            task1(argc, argv, a, b, n, func);
        } else {
            printf("Right integration limit is mandatory parameter.\n");
        }
    } else if(!strcmp(taskName, "2")) {
        task2(argc, argv, n, commMode);
    } else { // if invalid task number provided
        printf("Task %s doesn't exist. See help (-h) for available tasks.\n", taskName);
    }

    // end of program with exit code 0
    return 0;
}