#include <stdio.h>                          // import of the definitions of the C IO library
#include <string.h>                         // import of the definitions of the string operations
#include <unistd.h>                         // standard unix io library definitions and declarations
#include <errno.h>                          // system error numbers
#include <stdlib.h>                         // for random()
#include <time.h>                           // to seed random generator
#include <math.h>                           // for math functions
#include "mpi.h"                            // import of the MPI definitions

#include "util.h"                           // include own header file

void applyFilter(unsigned char input[][1280 + 4], unsigned char output[][1280], int filterType, int rowsPerProcessCount)
{
    double blurFilter[5][5] = {
        {0./37, 0./37, 1./37, 0./37, 0./37},
        {0./37, 2./37, 4./37, 2./37, 0./37},
        {1./37, 4./37, 9./37, 4./37, 1./37},
        {0./37, 2./37, 4./37, 2./37, 0./37},
        {0./37, 0./37, 1./37, 0./37, 0./37}
    };
    double sharpenFilter[5][5] = {
        {0, 0, 0, 0, 0},
        {0, 0, -1, 0, 0},
        {0, -1, 5, -1, 0},
        {0, 0, -1, 0, 0},
        {0, 0, 0, 0, 0}
    };
    double reliefFilter[5][5] = {
        {0, 0, 0, 0, 0},
        {0, -2, -1, 0, 0},
        {0, -1, 1, 1, 0},
        {0, 0, 1, 2, 0},
        {0, 0, 0, 0, 0}
    };
    double edgeFilter[5][5] = {
        {0./4, 0./4, 0./4, 0./4, 0./4},
        {0./4, 1./4, 2./4, 1./4, 0./4},
        {0./4, 2./4, -12./4, 2./4, 0./4},
        {0./4, 1./4, 2./4, 1./4, 0./4},
        {0./4, 0./4, 0./4, 0./4, 0./4},
    };

    switch(filterType){
        case 1: // blur
            getFilterResult(input, output, blurFilter, rowsPerProcessCount);
            break;
        case 2: // sharpen
            getFilterResult(input, output, sharpenFilter, rowsPerProcessCount);
            break;
        case 3: // relief
            getFilterResult(input, output, reliefFilter, rowsPerProcessCount);
            break;
        case 4: // edge
            getFilterResult(input, output, edgeFilter, rowsPerProcessCount);
            break;
        default:
            exit(1);
    }
}

void getFilterResult(unsigned char input[][1280 + 4], unsigned char output[][1280], double filter[][5], int rowsPerProcessCount)
{
    int x, y;   // pixel iterators

    for(y = 2; y < rowsPerProcessCount + 2; y++){
        for(x = 2; x < 1280 + 2; x++){
            output[y - 2][x - 2] = getSinglePixelFilterResult(input, filter, y, x);
        }
    }
}

unsigned char getSinglePixelFilterResult(unsigned char input[][1280 + 4], double filter[][5], int y, int x)
{
    int v, u,   // filter iterators
        k = 2;  // just for fun

    double result = 0.0;

    for(v = 0; v <= 2 * k; v++){
        for(u = 0; u <= 2 * k; u++){
            result += filter[v][u] * (double)(input[y + v - k][x + u - k]);
        }
    }

    if(result < 0.0){
        result = 0;
    } else if(result > 255.0){
        result = 255;
    }
    return (unsigned char)result;
}

void addLeftRightPadding(unsigned char input[][1280], unsigned char paddedInput[][1280 + 4], int rowsPerProcessCountPadded)
{
    int x, y;   // pixel iterators

    for(y = 0; y < rowsPerProcessCountPadded; y++){
        // add left padding
        for(x = 0; x < 2; x++){
            paddedInput[y][x] = 0;
        }

        // copy original image rows to the padded array
        for(x = 2; x < 1280 + 2; x++){
            paddedInput[y][x] = input[y][x - 2];
        }

        // add right padding
        for(x = 1280 + 2; x < 1280 + 4; x++){
            paddedInput[y][x] = 0;
        }
    }
}