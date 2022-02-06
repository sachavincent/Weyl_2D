
#include "utils.h"

int DEBUG = 0;

Image TransformMatrixToImage(Matrix *m)
{
    int height = MatNbRow(*m);
    int width = MatNbCol(*m);
    Image im = ImAlloc(GrayLevel, height, width);
    int **mat = MatGetInt(*m);
    unsigned char **mi = ImGetI(im);
#pragma omp for
    for (int x = 0; x < height; x++)
    {
        for (int y = 0; y < width; y++)
        {
            mi[x][y] = mat[x][y];
        }
    }
#pragma omp barrier

    if (DEBUG)
        printf("Converted Matrix Successfully\n");

    return im;
}

Matrix TransformImageToMatrix(Image *i)
{
    int height = ImNbRow(*i);
    int width = ImNbCol(*i);
    Matrix m = MatAlloc(Int, height, width);
    unsigned char **mat = ImGetI(*i);
    int **mi = MatGetInt(m);
#pragma omp for
    for (int x = 0; x < height; x++)
    {
        for (int y = 0; y < width; y++)
        {
            mi[x][y] = mat[x][y];
        }
    }
#pragma omp barrier

    if (DEBUG)
        printf("Converted Image Successfully\n");

    return m;
}
Matrix ReadMatrixFromPGM(char *imageName)
{
    Image image = ImRead(imageName);
    if (image == NULL)
        return NULL;

    Matrix mat = TransformImageToMatrix(&image);
    ImFree(&image);

    return mat;
}

Matrix GetIntegralImage(Matrix m)
{
    int height = MatNbRow(m);
    int width = MatNbCol(m);
    int **i = MatGetInt(m);

    Matrix mm = MatAlloc(Int, height, width);
    int **ii = MatGetInt(mm);

    ii[0][0] = i[0][0];
#pragma omp for
    for (int x = 1; x < height; x++)
    {
        ii[x][0] = ii[x - 1][0] + i[x][0];
    }

#pragma omp barrier
#pragma omp for
    for (int y = 1; y < width; y++)
    {
        int s = i[0][y];
        ii[0][y] = ii[0][y - 1] + s;
        for (int x = 1; x < height; x++)
        {
            s += i[x][y];
            ii[x][y] = ii[x][y - 1] + s;
        }
    }
#pragma omp barrier
    return mm;
}

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
Matrix SliceMatrix(Matrix M, int x, int y, int width, int height)
{
    assert(M != NULL);

    assert(width > 0 && height > 0);
    int heightM = MatNbRow(M);
    int widthM = MatNbCol(M);

    assert(widthM > 0 && heightM > 0);
    assert(x >= 0 && y >= 0 && x < widthM && y < heightM);

    int **values = MatGetInt(M);
    Matrix sliced = MatAlloc(Int, height, width);
    int **slicedValues = MatGetInt(sliced);

#pragma omp for
    for (int i = y; i < y + height; i++)
    {
        for (int j = x; j < x + width; j++)
        {
            slicedValues[i - y][j - x] = values[i][j];
        }
    }
#pragma omp barrier

    if (DEBUG)
        printf("Successfully sliced matrix!\n");
    return sliced;
}
/**
 * @brief Get the Sliced Difference Matrix object
 *
 * @param A whole matrix
 * @param x top left x coordinate
 * @param y top left y coordinate
 * @param B smaller matrix
 *
 * @return difference between A and B at (x, y)
 */
Matrix GetSlicedDifferenceMatrix(Matrix A, int x, int y, Matrix B)
{
    int h = MatNbRow(B);
    int w = MatNbCol(B);
    int **a = MatGetInt(A);
    int **b = MatGetInt(B);
    Matrix sliced = MatAlloc(Int, h, w);
    int **m = MatGetInt(sliced);

#pragma omp for
    for (int i = y; i < y + h; i++)
    {
        for (int j = x; j < x + w; j++)
        {
            m[i - y][j - x] = a[i][j] - b[i - y][j - x];
        }
    }

#pragma omp barrier
    return sliced;
}

int WeylA(int **pi, int x, int y, int W, int H)
{
    return pi[H][W] + pi[x][y] - pi[H][y] * pi[H][y] - pi[x][W];
}

int *WeylFourSquare(Matrix m, int x, int y)
{
    int *allPi = (int *)malloc(sizeof(int) * 4);
    int **pi = MatGetInt(m);

    int H = MatNbRow(m) - 1;
    int W = MatNbCol(m) - 1;

    int pi1, pi2, pi3, pi4;

    pi1 = pi[x][y];
    if (x != H)
    {
        pi2 = pi[H][y] - pi1;
    }
    else
    {
        pi2 = pi[H][y];
    }

    if (y != W)
    {
        pi3 = pi[x][W] - pi1;
    }
    else
    {
        pi3 = pi[x][W];
    }
    if (x != H)
    {
        if (y != W)
        {
            pi4 = WeylA(pi, x, y, W, H);
        }
        else
        {
            pi4 = pi[H][y] - pi1;
        }
    }
    else
    {
        if (y != W)
        {
            pi4 = pi[x][W] - pi1;
        }
        else
        {
            pi4 = pi1;
        }
    }
    // voir (4) et (5) papier On a fast implementation ..
    allPi[0] = pi1;
    allPi[1] = pi2;
    allPi[2] = pi3;
    allPi[3] = pi4;
    return allPi;
}

int WeylNorm(Matrix m)
{
    Matrix integralImage = GetIntegralImage(m);
    int H = MatNbRow(m);
    int W = MatNbCol(m);

    int minPi1 = INT_MAX, minPi2 = INT_MAX, minPi3 = INT_MAX, minPi4 = INT_MAX;
    int maxPi1 = INT_MIN, maxPi2 = INT_MIN, maxPi3 = INT_MIN, maxPi4 = INT_MIN;

#pragma omp for
    for (int x = 0; x < H; x++)
    {
        for (int y = 0; y < W; y++)
        {
            int *allPi = WeylFourSquare(integralImage, x, y);
            minPi1 = min(minPi1, allPi[0]);
            minPi2 = min(minPi2, allPi[1]);
            minPi3 = min(minPi3, allPi[2]);
            minPi4 = min(minPi4, allPi[3]);

            maxPi1 = max(maxPi1, allPi[0]);
            maxPi2 = max(maxPi2, allPi[1]);
            maxPi3 = max(maxPi3, allPi[2]);
            maxPi4 = max(maxPi4, allPi[3]);

            free(allPi);
        }
    }
#pragma omp barrier
    int pi1, pi2, pi3, pi4;
    pi1 = maxPi1 - minPi1;
    pi2 = maxPi2 - minPi2;
    pi3 = maxPi3 - minPi3;
    pi4 = maxPi4 - minPi4;
    MatFree(&integralImage);
    return max(max(pi1, pi2), max(pi4, pi3));
}

int *WeylMinMaxOpti(int pi, int c)
{
    int *minMax = (int *)malloc(sizeof(int) * 2);
    minMax[0] = max(c, c - pi);
    minMax[1] = min(c, c - pi);
    return minMax;
}

Matrix GetIntegralImageOpti(Matrix m, int *maxPerRow, int *minPerRow)
{
    int height = MatNbRow(m);
    int width = MatNbCol(m);
    int **i = MatGetInt(m);

    Matrix mm = MatAlloc(Int, height, width);

    int **ii = MatGetInt(mm);

    ii[0][0] = i[0][0];
    maxPerRow[0] = ii[0][0];
    minPerRow[0] = ii[0][0];
#pragma omp for
    for (int x = 1; x < height; x++)
    {
        ii[x][0] = ii[x - 1][0] + i[x][0];
        maxPerRow[x] = ii[x][0];
        minPerRow[x] = ii[x][0];
    }

#pragma omp barrier
#pragma omp for
    for (int y = 1; y < width; y++)
    {
        int s = i[0][y];
        ii[0][y] = ii[0][y - 1] + s;
        maxPerRow[0] = max(ii[0][y], maxPerRow[0]);
        minPerRow[0] = min(ii[0][y], minPerRow[0]);
        for (int x = 1; x < height; x++)
        {
            s += i[x][y];
            ii[x][y] = ii[x][y - 1] + s;
            maxPerRow[x] = max(ii[x][y], maxPerRow[x]);
            minPerRow[x] = min(ii[x][y], minPerRow[x]);
        }
    }
#pragma omp barrier
    return mm;
}

int WeylNormOpti(Matrix m)
{
    int H = MatNbRow(m);
    int W = MatNbCol(m);
    int *minPerRow = (int *)malloc(sizeof(int) * H);
    int *maxPerRow = (int *)malloc(sizeof(int) * H);
    Matrix integralImage = GetIntegralImageOpti(m, maxPerRow, minPerRow);

    int **ii = MatGetInt(integralImage);
    int minPi = INT_MAX;
    int maxPi = INT_MIN;
#pragma omp for
    for (int x = 0; x < H; x++)
    {
        int miniRow = minPerRow[x];
        int maxiRow = maxPerRow[x];
        // printf("min : %d max : %d\n", miniRow, maxiRow);
        int c = ii[x][W];
        maxPi = max(maxPi, max(c, c - miniRow));
        minPi = min(minPi, min(c, c - maxiRow));
        /*
        for (int y = 0; y < W; y++)
        {
            // int *allPi = WeylMinMaxOpti(ii[x][y], c);m
            maxPi = max(maxPi, max(c, c - ii[x][y]));
            minPi = min(minPi, min(c, c - ii[x][y]));
            // free(allPi);
        }
        */
    }
#pragma omp barrier
    MatFree(&integralImage);
    free(minPerRow);
    free(maxPerRow);
    // printf("Test : maxPi = %d\n", maxPi);
    return maxPi - minPi;
}

/**
 * @brief Get neighborhood pixels in Matrix centered on (xC, yC) of size neigh*neigh
 *
 * @param M Matrix of pixels
 * @param xC center of neighborhood
 * @param yC center of neighborhood
 * @param neigh size of neighborhood (must be odd)
 * @return Matrix containing the neighborhood
 */
Matrix GetNeighborhood(Matrix M, int xC, int yC, int neigh)
{
    assert(neigh % 2 == 1);
    assert(neigh > 0);

    int **m = MatGetInt(M);

    int halfNeigh = neigh / 2;

    // printf("%d %d\n", m[xC][yC], halfNeigh);
    Matrix res = MatAlloc(Int, neigh, neigh);
    int **mRes = MatGetInt(res);
#pragma omp for
    for (int x = xC - halfNeigh; x <= xC + halfNeigh; x++)
    {
        for (int y = yC - halfNeigh; y <= yC + halfNeigh; y++)
        {
            int pixel = m[x][y];

            mRes[x - (xC - halfNeigh)][y - (yC - halfNeigh)] = pixel;
        }
    }

#pragma omp barrier
    return res;
}