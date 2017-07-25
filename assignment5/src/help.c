#include "help.h"                       // own header file
#include <stdio.h>                      // import of the definitions of the C IO library
#include "mpi.h"                        // import of the MPI definitions

void showHelp(int assignmentNumber, int argc, char* argv[])
{
    int myRank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    if(myRank == 0){
        printf("#########################################################################################\n");
        printf("# Assignment %d                                                                          #\n", assignmentNumber);
        printf("#########################################################################################\n");
        printf("# To compile the program execute in root directory:                                     #\n");
        printf("#  -> make clean && make                                                                #\n");
        printf("#                                                                                       #\n");
        printf("# To run the program on all available hosts execute in root directory:                  #\n");
        printf("#  -> mpiexec ./assignment [arguments]                                                  #\n");
        printf("#                                                                                       #\n");
        printf("# You can specify the amount of hosts by using -np argument:                            #\n");
        printf("#  -> mpiexec -np 4 ./assignment [arguments]                                            #\n");
        printf("#                                                                                       #\n");
        printf("#########################################################################################\n");
        printf("#                                                                                       #\n");
        printf("# Arguments:                                                                            #\n");
        printf("#                                                                                       #\n");
        printf("# -h                   view this help                                                   #\n");
        printf("# -t [integer]         task number to execute                                           #\n");
        printf("#                      valid task numbers are 1 and 2                                   #\n");
        printf("#                                                                                       #\n");
        printf("# Available parameters for task 1                                                       #\n");
        printf("#                                                                                       #\n");
        printf("# -i [string]          input image to which filter will be applied                      #\n");
        printf("#                      if not specified, default value is ./examples/ffm_1280x960.gray  #\n");
        printf("# -f [integer]         filter type; valid filter types are:                             #\n");
        printf("#                      1 = blur (default), 2 = sharpen, 3 = relief, 4 = edge            #\n");
        printf("# -m [integer]         filter strength = how often to apply filter                      #\n");
        printf("#                      if not specified, default value is 1                             #\n");
        printf("#                                                                                       #\n");
        printf("# Available parameters for task 2                                                       #\n");
        printf("#                                                                                       #\n");
        printf("# None. Task 2 is interactive.                                                          #\n");
        printf("#                                                                                       #\n");
        printf("#########################################################################################\n");
    }
    // finalizing MPI interface
    MPI_Finalize();
}