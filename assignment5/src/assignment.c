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
#include "task2.h"                      // to run Task 1

// constants & default values
#define ASSIGNMENT_NR       5
#define INPUT_IMG_FILE      ((unsigned char *)"./examples/ffm_1280x960.gray")    // default input image path
int main(int argc, char* argv[])
{
    int opt = 0,                            // cli parameter name
        taskNumber = -1,                    // cli parameter for task number, default 0 (not selected)
        filterStrength = 1,                 // cli parameter for filter strength, default 1
        filterType = 1;                     // cli parameter for filter type, default blur

    char *inputImgFilePath = INPUT_IMG_FILE;

    // iterate through all cli parameters
    while ((opt = getopt(argc, argv, "i:m:f:t:h")) != -1) {
        switch(opt) {
            case 'h': // if -h parameter given => show help
                showHelp(ASSIGNMENT_NR, argc, argv);
                exit(0);
            case 't': // if -t parameter given => set task number to execute
                taskNumber = strtoumax(optarg, NULL, 10);
                break;
            case 'i': // if -i parameter given => set path to input image file
                inputImgFilePath = optarg;
                break;
            case 'f': // if -f parameter given => set which filter to apply
                filterType = strtoumax(optarg, NULL, 10);
                break;
            case 'm': // if -m parameter given => set how often to apply filter
                filterStrength = strtoumax(optarg, NULL, 10);
                break;
            default: // if any other cli parameters provided
                // exit program with error status
                exit(1);
        }
    }

    // decide which task to run based on -t parameter
        if(taskNumber == -1){ // if no parameter provided - just show help
            showHelp(ASSIGNMENT_NR, argc, argv);
            exit(0);
        } else if(taskNumber == 1) {
            task1(argc, argv, inputImgFilePath, filterType, filterStrength);
        } else if(taskNumber == 2) {
            task2(argc, argv);
        } else { // if invalid task number provided
            printf("Task %s doesn't exist. See help (-h) for available tasks.\n", taskNumber);
        }

    // end of program with exit code 0
    return 0;
}