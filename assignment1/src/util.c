#include <stdio.h> 	    // import of the definitions of the C IO library
#include <string.h>     // import of the definitions of the string operations
#include <unistd.h>	    // standard unix io library definitions and declarations
#include <errno.h>	    // system error numbers
#include <stdlib.h> 	// for random()
#include <time.h>		// to seed random generator

#include "util.h" 		// include own header file

void fillWithRandomInt(int* array, int size, int maxRandom) 
{
		srandom( (unsigned) time(NULL) );
	int i;
	for(i = 0; i < size; i++) {
		// assign random int to each array position
		array[i] = generateRandomInt(maxRandom);
	} 
}

int getMin(int* array, int size)
{
    int i;
    int min = array[0];

 	for (i = 1; i < size; i++) {
        if (array[i] < min) {
           min = array[i];
	   	}
    }
	return min;
}

int getMax(int* array, int size)
{
    int i;
    int max = array[0];

 	for (i = 1; i < size; i++) {
        if (array[i] > max) {
           max = array[i];
	   	}
    }
	return max;
}

int getSum(int* array, int size)
{
    int i;
    int sum = 0;

 	for (i = 0; i < size; i++) {
		sum += array[i];
	}
	return sum;
}

int generateRandomInt(int max){
	return random() % max;
}