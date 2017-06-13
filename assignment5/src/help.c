#include "help.h"                       // own header file
#include <stdio.h>                      // import of the definitions of the C IO library
#include "mpi.h"                        // import of the MPI definitions

void showHelp(int assignmentNumber, int argc, char* argv[])
{
    int myRank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    if(myRank == 0){
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
        printf("# -h                   view this help                                     #\n");
        printf("# -m [integer]         amount of random integers per process              #\n");
        printf("#                      if not specified, default value is 1               #\n");
        printf("###########################################################################\n");
    }
    // finalizing MPI interface
    MPI_Finalize();
}