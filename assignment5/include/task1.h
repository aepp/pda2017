#ifndef TASK1_H_INCLUDED
#define TASK1_H_INCLUDED
/* ^^ these are the include guards */

/* Prototypes for the functions */

/**
 * Task 1 of Assignment 5
 *
 * Parameters:
 *
 * argc                     number of cli arguments
 * argv                     cli arguments
 * inputImgFilePath         path to input image file
 * filterType               filter type to apply (see help for details)
 * filterStrength           how often to apply desired filter
 */
void task1(int argc, char* argv[], char* inputImgFilePath, int filterType, int filterStrength);

#endif