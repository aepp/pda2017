#include "help.h"                       // own header file
#include <stdio.h>                      // import of the definitions of the C IO library
#include "mpi.h"                        // import of the MPI definitions

void showHelp(int assignmentNumber, double epsilon, int argc, char* argv[])
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
        printf("# -t [integer]         run task specified by [integer] from Assignment %d #\n", assignmentNumber);
        printf("#                      valid task numbers are 1 and 2                     #\n");
        printf("# -h                   view this help                                     #\n");
        printf("#                      if not set the default value is 100                #\n");
        printf("# -m [string]          path to matrix file                                #\n");
        printf("#                      if not specified, set to default value             #\n");
        printf("# -v [string]          path to vector file                                #\n");
        printf("#                      if not specified, set to default value             #\n");
        printf("# -s [integer]         size of vector b                                   #\n");
        printf("#                      if not specified, the default value is 8           #\n");
        printf("# -e [double]          epsilon value, denotes convergence value           #\n");
        printf("#                      for Jacobi method                                  #\n");
        printf("#                      if not set the default value is %lf           #\n", epsilon);
        printf("###########################################################################\n");
    }
    // finalizing MPI interface
    MPI_Finalize();
}