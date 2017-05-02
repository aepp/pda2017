/************************************************************************/
/* Program: main.c		                                                */ 
/* Author: Aleksandr Epp <aleksandr.epp@gmail.com>						*/
/* matriclenumber: 6002853                                              */
/* Assignment : 1                                                       */	
/* Parameters: -h, --help, -m                                           */
/* Environment variables: no                                            */
/*                                                                      */
/* Description:                                                         */
/*                                                                      */
/* assignment1 															*/
/*                                                                      */
/************************************************************************/ 

#include <stdio.h> 	    // import of the definitions of the C IO library
#include <string.h>     // import of the definitions of the string operations
#include <unistd.h>	    // standard unix io library definitions and declarations
#include <errno.h>	    // system error numbers

#include "help.h"
#include "task1.h"
#include "task2a.h"
#include "task2b.h"

int main(int argc, char* argv[ ]) 
{ 
	int assignmentNumber = 1,
		maxRandom = 100; // maximum random integer that appears in the array
	
	if (argc == 2) {
		if (!strcmp(argv[1], "--help") || !strcmp(argv[1], "-h")) {
			showHelp(assignmentNumber);
		} else if(!strcmp(argv[1], "--task1")) {
			task1(&argc, &argv);
		} else if(!strcmp(argv[1], "--task2a")) {
			task2a(&argc, &argv, maxRandom);
		} else if(!strcmp(argv[1], "--task2b")) {
			task2b(&argc, &argv);
		} 
	} else {
		showHelp(assignmentNumber);
	}
	// end of progam with exit code 0 
	return 0;
}