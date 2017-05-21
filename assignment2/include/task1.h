#ifndef TASK1_H_INCLUDED
#define TASK1_H_INCLUDED
/* ^^ these are the include guards */

/* Prototypes for the functions */

/**
 * Task 1 of Assignment 2
 *
 * Parameters:
 *
 * argc         number of cli arguments
 * argv         cli arguments
 * a            left integration limit (default = 1)
 * b            right integration limit
 * n            number of intervals  (default = 100)
 * funcNumber   function to integrate (default = 1)
 */
void task1(int argc, char* argv[], double a, double b, double n, double funcNumber);

#endif