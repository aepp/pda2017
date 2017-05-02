/************************************************************************/
/* Author: Aleksandr Epp <aleksandr.epp@gmail.com>                      */
/* matriclenumber: 6002853                                              */
/* Assignment : 1                                                       */
/* Task : 1                                                             */
/*                                                                      */
/* Description:                                                         */
/*                                                                      */
/* For each CPU, where the program is running, the rank and the host    */
/* name are sent to the process with rank 0. Process with rank 0 then   */
/* prints received messages.                                            */
/*                                                                      */
/************************************************************************/ 

#include "mpi.h"        // import of the MPI definitions
#include <stdio.h>      // import of the definitions of the C IO library
#include <string.h>     // import of the definitions of the string operations
#include <unistd.h>     // standard unix io library definitions and declarations
#include <errno.h>      // system error numbers

#include "task1.h"      // include own header file

void task1(int* argc, char** argv[ ]) {
    // length of name
    int namelen;
    // rank of the process
    int my_rank;
    // hostname
    char *c, proc_name[MPI_MAX_PROCESSOR_NAME+1];
    
    // rank of the process sent the message
    int theirRank;
    // name of the process sent the message
    char theirProcName[MPI_MAX_PROCESSOR_NAME+1];
    
    // default tag
    int MY_TAG = 0;
    int TAG_PROC_NAME = 1;
    int TAG_PROC_NUMB = 2;
    // number of proccesses envolved
    int numProc; 
    // iterator variable
    int i;
    
    // initializing of MPI-Interface
    MPI_Init(argc, argv);
    //get your rank
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // get the number of processes running
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);

    // initialize string with NULL-characters
    memset (proc_name, 0, MPI_MAX_PROCESSOR_NAME+1 );
                                    
    // finding out own computer name
    MPI_Get_processor_name(proc_name, &namelen);
    
    // separate the first part of hostname
    if ( (c=strchr(proc_name,'.'))  !=  NULL) *c = '\0';     
    
    //printf("%03d: process running on %s!\n", my_rank, proc_name );
                                    // print rank and hostname to display
    
    if(my_rank == 0){ // if I'm process with rank 0
        // recieve messages from all other processes
        for(i = 0; i < numProc - 1; i++){
            // recieve the process rank
            MPI_Recv(&theirRank, 1, MPI_INT, MPI_ANY_SOURCE, TAG_PROC_NUMB, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // recieve the process name
            MPI_Recv(&theirProcName, MPI_MAX_PROCESSOR_NAME+1, MPI_CHAR, MPI_ANY_SOURCE, TAG_PROC_NAME, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // print the recieved message to console
            printf("Msg recieved from %s: %d!\n", theirProcName, theirRank);
        }
    } else { // if I'm process with rank > 0, send message to process with rank 0
        // send my rank
        MPI_Send(&my_rank, 1, MPI_INT, 0, TAG_PROC_NUMB, MPI_COMM_WORLD);
        // send my name
        MPI_Send(&proc_name, namelen+1, MPI_CHAR, 0, TAG_PROC_NAME, MPI_COMM_WORLD);
    }
    // finalizing MPI interface 
    MPI_Finalize();                    
}