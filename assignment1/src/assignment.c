/************************************************************************/
/* Program: main.c                                                      */
/* Author: Aleksandr Epp <aleksandr.epp@gmail.com>                      */
/* matriclenumber: 6002853                                              */
/* Assignment : 1                                                       */
/* Parameters: -h, -t, -m                                               */
/* Environment variables: no                                            */
/*                                                                      */
/* Description:                                                         */
/*                                                                      */
/* assignment1                                                          */
/*                                                                      */
/************************************************************************/ 

#include <stdio.h>      // import of the definitions of the C IO library
#include <stdlib.h>     // to use exit() function
#include <string.h>     // import of the definitions of the string operations
#include <unistd.h>     // standard unix io library definitions and declarations
#include <errno.h>      // system error numbers
#include <inttypes.h>   // for strtoumax() string to int coversion

#include "help.h"
#include "task1.h"
#include "task2a.h"
#include "task2b.h"

#define DEFAULT_ARRAY_SIZE 100 // default number of elements in random int array
#define MAX_RANDOM 100 // maximum random integer that appears in the array
#define ASSIGNMENT_NR 1

int main(int argc, char* argv[ ]) 
{ 
    int arraySize = DEFAULT_ARRAY_SIZE, // number of elements in random int array
        i, // iteration variable
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
//                printf("Task will run with %d elements\n", arraySize);
                break;
            case 't': // if -t parameter given => save task nr. to decide which task to execute
                taskName = optarg;
//                printf("Execute task %s\n", taskName);
                break;
            default: // if any other cli parameters provided show
                //printf("\nInvalid parameter or parameter value missing\n");
                exit(1);
        }
    }

    // if no parameter provided - just show help
    if(taskName == NULL){
        showHelp(ASSIGNMENT_NR);
        exit(0);
    } else if(!strcmp(taskName, "1")) {
        task1(&argc, &argv);
    } else if(!strcmp(taskName, "2a")) {
        task2a(argc, argv, MAX_RANDOM, arraySize);
    } else if(!strcmp(taskName, "2b")) {
        task2b(argc, argv, MAX_RANDOM, arraySize);
    } else {
        printf("Task %s doesn't exist. See help (-h) for available tasks.\n", taskName);
    }

    // end of program with exit code 0
    return 0;
}