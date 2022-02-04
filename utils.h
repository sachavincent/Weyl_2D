
#ifndef __utils_h_
#define __utils_h_

#include <omp.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <getopt.h>
#include <sys\timeb.h>
#include <errno.h>
#include <limits.h>
#include "limace.h"

#define min(x, M) ((x) >= M ? M : x)
#define max(x, M) ((x) < M ? M : x)

extern int DEBUG;

extern Image TransformMatrixToImage(Matrix *m);

extern Matrix TransformImageToMatrix(Image *i);

extern Matrix ReadMatrixFromPGM(char *imageName);

extern Matrix GetIntegralImage(Matrix m);

extern Matrix GetSlicedDifferenceMatrix(Matrix A, int x, int y, Matrix B);

extern int WeylA(int **pi, int x, int y, int W, int H);

extern int *WeylFourSquare(Matrix m, int x, int y);

extern int WeylNorm(Matrix m);

extern int *WeylMinMaxOpti(int pi, int c);

extern Matrix GetIntegralImageOpti(Matrix m, int *maxPerRow, int *minPerRow);

extern int WeylNormOpti(Matrix m);

/**
 * @brief Get neighborhood pixels in Matrix centered on (xC, yC) of size neigh*neigh
 *
 * @param M Matrix of pixels
 * @param xC center of neighborhood
 * @param yC center of neighborhood
 * @param neigh size of neighborhood (must be odd)
 * @return Matrix containing the neighborhood
 */
extern Matrix GetNeighborhood(Matrix M, int xC, int yC, int neigh);

/**
 * @brief Slice matrix
 *
 * @param M original matrix
 * @param x top left x coordinate
 * @param y top left y coordinate
 * @param width width of slice
 * @param height height of slice
 *
 * @return result
 */
Matrix SliceMatrix(Matrix M, int x, int y, int width, int height);

#endif