#ifndef UTIL_H_INCLUDED
#define UTIL_H_INCLUDED
/* ^^ these are the include guards */

/* Prototypes for the functions */

/**
 * apply selected filter
 *
 * input        part of input image
 * output       part of filtered output result
 * filterType   filter selector
 */
void applyFilter(unsigned char input[][1280 + 4], unsigned char output[][1280], int filterType, int rowsPerProcessCount);

/**
 * calculate filtered values
 *
 * input        part of input image
 * output       part of filtered output result
 * filter       selected filter (matrix)
 */
void getFilterResult(unsigned char input[][1280 + 4], unsigned char output[][1280], double filter[][5], int rowsPerProcessCount);

/**
 * calculate filter result for a single pixel
 *
 * input        part of input image
 * filter       selected filter (matrix)
 * y            y position of current pixel
 * x            x position of current pixel
 */
unsigned char getSinglePixelFilterResult(unsigned char input[][1280 + 4], double filter[][5], int y, int x);

/**
 * add left and right padding of 2 pixels to an image part
 *
 * input                        part of input image (padded at top and bottom)
 * paddedInput                  part of input image after left and right padding added
 * rowsPerProcessCountPadded    rows count in the completely padded image
 */
void addLeftRightPadding(unsigned char input[][1280], unsigned char paddedInput[][1280 + 4], int rowsPerProcessCountPadded);

#endif