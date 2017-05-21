#include <stdio.h>          // import of the definitions of the C IO library
#include <string.h>         // import of the definitions of the string operations
#include <unistd.h>         // standard unix io library definitions and declarations
#include <errno.h>          // system error numbers
#include <stdlib.h>         // for random()
#include <time.h>           // to seed random generator
#include <math.h>           // for math functions

#include "util.h"           // include own header file

double f1(double x){
    return 1 / (sqrt(2 * x + 1));
}

double f2(double x){
    return log(x);
}

double trapezoidalRuleF1(double a, double b, double s, int n){
    int k;
    double tempSum = 0;

    for(k = 1; k < n; k++) {
        tempSum += f1(a + k * s);
    }
    return (s/2) * (f1(a) + f1(b) + 2 * tempSum);
}

double trapezoidalRuleF2(double a, double b, double s, int n){
    int k;
    double tempSum = 0;

    for(k = 1; k < n; k++) {
        tempSum += f2(a + k * s);
    }
    return (s/2) * (f2(a) + f2(b) + 2 * tempSum);
}

void fillWithRandomInt(int* array, int size, int maxRandom, int rank)
{
    srandom((unsigned)time(NULL) * rank);
    int i;
    for(i = 0; i < size; i++) {
        // assign random int to each array position
        array[i] = generateRandomInt(maxRandom);
    }
}

int generateRandomInt(int max){
    return random() % max;
}
