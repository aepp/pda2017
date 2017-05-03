#include "help.h"       // own header file
#include <stdio.h>      // import of the definitions of the C IO library

void showHelp(int assignmentNumber)
{
    printf("#########################################################################\n");
    printf("# Assignment %d                                                          #\n", assignmentNumber);
    printf("#########################################################################\n");
    printf("# To compile the program execute in root directory:                     #\n");
    printf("#  -> make clean && make                                                #\n");
    printf("#                                                                       #\n");
    printf("# To run the program on all available hosts execute in root directory:  #\n");
    printf("#  -> mpiexec ./assignment [arguments]                                  #\n");
    printf("#                                                                       #\n");
    printf("# You can specify the amount of hosts by using -np argument:            #\n");
    printf("#  -> mpiexec -np 4 ./assignment [parameter]                            #\n");
    printf("#                                                                       #\n");
    printf("#########################################################################\n");
    printf("# Arguments:                                                            #\n");
    printf("#                                                                       #\n");
    printf("# -t (task number)      run task (task number) from Assignment %d        #\n", assignmentNumber);
    printf("#                       available tasks are 1, 2a and 2b                #\n");
    printf("# -h                    view this help                                  #\n");
    printf("# -m [array size]       number of elements in random array              #\n");
    printf("#                       (if not set use default 100)                    #\n");
    printf("#########################################################################\n");
}