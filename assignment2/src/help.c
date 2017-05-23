#include "help.h"       // own header file
#include <stdio.h>      // import of the definitions of the C IO library

void showHelp(int assignmentNumber)
{
    printf("###########################################################################\n");
    printf("# Assignment %d                                                            #\n", assignmentNumber);
    printf("###########################################################################\n");
    printf("# To compile the program execute in root directory:                       #\n");
    printf("#  -> make clean && make                                                  #\n");
    printf("#                                                                         #\n");
    printf("# To run the program on all available hosts execute in root directory:    #\n");
    printf("#  -> mpiexec ./assignment [arguments]                                    #\n");
    printf("#                                                                         #\n");
    printf("# You can specify the amount of hosts by using -np argument:              #\n");
    printf("#  -> mpiexec -np 4 ./assignment [arguments]                              #\n");
    printf("#                                                                         #\n");
    printf("###########################################################################\n");
    printf("# Arguments:                                                              #\n");
    printf("#                                                                         #\n");
    printf("# -t [task number]      run task [task number] from Assignment %d          #\n", assignmentNumber);
    printf("#                       valid task numbers are 1 and 2                    #\n");
    printf("# -h                    view this help                                    #\n");
    printf("# -n [integer]          (task 1) number of intervals                      #\n");
    printf("#                       (task 2) amount of random numbers                 #\n");
    printf("#                       if not set the default value is 100               #\n");
    printf("# -a [left limit]       (task 1 only) left integration limit              #\n");
    printf("#                       if not set the default value is 0                 #\n");
    printf("# -b [right limit]      (task 1 only) right integration limit             #\n");
    printf("# -f [function number]  (task 1 only) select which function to integrate  #\n");
    printf("#                       valid functions are 1 (=sqrt) and 2 (=log)        #\n");
    printf("#                       if not set the default value is 1                 #\n");
    printf("# -m [mode]             (task 2 only) communication mode                  #\n");
    printf("#                       valid modes are 1 (=task2a) and 2 (=task2b)       #\n");
    printf("#                       if not set the default value is 1                 #\n");
    printf("###########################################################################\n");
}