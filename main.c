#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "limace.h"
#include <assert.h>
#include <omp.h>
#include <sys\timeb.h>

#define min(x, M) ((x) >= M ? M : x)
#define max(x, M) ((x) < M ? M : x)

int DEBUG = 1;

int distance(int a, int b)
{
	return abs(b - a);
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

Matrix DifferenceMatrices(Matrix A, Matrix B)
{
	int h = MatNbRow(B);
	int w = MatNbCol(B);
	Matrix diff = MatAlloc(Int, h, w);
	int **a = MatGetInt(A);
	int **b = MatGetInt(B);
	int **m = MatGetInt(diff);
#pragma omp for
	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
		{
			m[i][j] = a[i][j] - b[i][j];
		}
	}
#pragma omp barrier
	return diff;
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
Matrix Stereo(Matrix Im1, Matrix Im2, int neigh, int seuilPerc, int opti)
{
	assert(seuilPerc >= 0 && seuilPerc <= 100);

	int width1 = MatNbCol(Im1);
	int height1 = MatNbRow(Im1);
	int width2 = MatNbCol(Im2);
	int height2 = MatNbRow(Im2);

	assert(width1 > 0 && height1 > 0);
	assert(width1 == width2 && height1 == height2);

	int halfNeigh = neigh / 2;

	Matrix disparityMap = MatAlloc(Int, height1, width1);
	int **mapValues = MatGetInt(disparityMap);
	Matrix norms = MatAlloc(Int, height1, width1);
	int **normsValues = MatGetInt(norms);

	int maxValue = INT_MIN;
	int maxNorm = INT_MIN;
	Matrix toGetSize = MatAlloc(Int, neigh, neigh);
	// Parcours de Im1

	for (int x = halfNeigh; x < height1 - halfNeigh; x++)
	{
		Matrix *listImagettes = (Matrix *)malloc(sizeof(toGetSize) * (width1 - halfNeigh * 2));
#pragma omp for
		for (int y = halfNeigh; y < width1 - halfNeigh; y++)
		{
			listImagettes[y - halfNeigh] = GetNeighborhood(Im2, x, y, neigh);
		}
#pragma omp barrier
#pragma omp for
		for (int y = halfNeigh; y < width1 - halfNeigh; y++)
		{
			int bestNorm = INT_MAX;
			Matrix imagette = GetNeighborhood(Im1, x, y, neigh);
			// Recherche dans Im2
			int bestY;
			int distance = width1 - halfNeigh - y;
			for (int y2 = y; y2 < width1 - halfNeigh; y2++)
			{
				Matrix diff = DifferenceMatrices(imagette, listImagettes[y2 - halfNeigh]);
				int norm;
				if (opti)
					norm = WeylNormOpti(diff);
				else
					norm = WeylNorm(diff);
					
				if (bestNorm > norm)
				{
					bestNorm = norm;
					bestY = y2;
				}

				MatFree(&diff);
			}
			MatFree(&imagette);
			normsValues[x][y] = bestNorm;
			maxNorm = max(maxNorm, bestNorm);
			mapValues[x][y] = (int)(255.0 * ((double)abs(bestY - y) / (double)distance));
			if (maxValue < mapValues[x][y])
				maxValue = mapValues[x][y];
		}

#pragma omp barrier

#pragma omp for
		for (int y = halfNeigh; y < width1 - halfNeigh; y++)
		{
			MatFree(&listImagettes[y - halfNeigh]);
		}

#pragma omp barrier
		free(listImagettes);
	}
	printf("Max norm: %d, perc:%d\n", maxNorm, seuilPerc);
	double seuil = ((double)(100 - (double)seuilPerc) / 100.0) * (double)maxNorm;
	printf("Seuil: %f\n", seuil);

#pragma omp for
	for (int x = 0; x < height1; x++)
	{
		for (int y = 0; y < width1; y++)
		{
			if (normsValues[x][y] > seuil)
				mapValues[x][y] = 0;
		}
	}

#pragma omp barrier

	return disparityMap;
}

Matrix GetSlicedDifferenceMatrix(Matrix A, int x, int y, Matrix B)
{
	int h = MatNbRow(B);
	int w = MatNbCol(B);
	int **m = MatGetInt(B);
	Matrix sliced = MatAlloc(Int, h, w);
	int **a = MatGetInt(A);
	int **b = MatGetInt(sliced);
#pragma omp for
	for (int i = x; i < x + h; i++)
	{
		//#pragma omp for
		for (int j = y; j < y + w; j++)
		{
			m[i - x][j - y] = a[i][j] - b[i - x][j - y];
		}
	}
#pragma omp barrier
	return sliced;
}

Matrix FindBInA(Matrix A, Matrix B)
{
	int heightA = MatNbRow(A);
	int widthA = MatNbCol(A);

	int heightB = MatNbRow(B);
	int widthB = MatNbCol(B);
	if (heightB > heightA || widthB > widthA)
	{
		printf("Error in findBInA, image B is larger than A\n");
		return MatAlloc(Int, 0, 0);
	}
	int normMin = INT_MAX;
	Matrix found = GetSlicedDifferenceMatrix(A, 0, 0, B);
	// int **foundValues = MatGetInt(found);
	Matrix disparityMatrix = MatAlloc(Int, heightA - heightB, widthA - widthB);
	int **dm = MatGetInt(disparityMatrix);
	// Matrix slicedA = GetSlicedMatrix(A, 0, 0, heightB, widthB);

#pragma omp for
	for (int x = 0; x < heightA - heightB; x++)
	{
		////#pragma omp for
		for (int y = 0; y < widthA - widthB; y++)
		{
			Matrix differenceMatrix = GetSlicedDifferenceMatrix(A, x, y, B);
			int newNorm = WeylNorm(differenceMatrix);
			dm[x][y] = newNorm;
			if (newNorm <= normMin)
			{
				normMin = newNorm;
				found = MatCopy(differenceMatrix);
			}
			MatFree(&differenceMatrix);
		}
	}
#pragma omp barrier

	MatWriteAsc(found, "foundMatrix.mx");
	MatWriteAsc(disparityMatrix, "disparityMatrix.mx");
	MatFree(&disparityMatrix);
	return found;
}
Matrix transformImageToMatrix(Image *i)
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

	printf("Converted Image Success\n");
	return m;
}

int main(int argc, char *argv[])
{
	if (argc != 3)
	{
		printf("Error parameter");
		return 1;
	}
	/*
		Image i = ImRead(argv[1]);
		if (i == NULL)
		{
			printf("Probleme lors de la lecture de l'image");
			return 2;
		}
		Matrix m = transformImageToMatrix(&i);
		//Matrix m = MatReadAsc(argv[1]);
		if (m == NULL)
		{
			printf("Probleme lors de la lecture de la matrice");
			return 2;
		}
		//Matrix M_IG = GetIntegralImage(m);
		//MatWriteAsc(M_IG, "test.mx");
		//MatWriteAsc(m,"imagetteCarre.mx");
		//int D = WeylNorm(m);
		//printf("Norme attendue=193917\nNorme obtenue=%d\n", D);

		//	Matrix integralImage = getIntegralImage(m);
		//	double dist = WeylDistance(integralImage);
		// Norme attendue = 193917

		Image realImage = ImRead("allCoins.pgm");
		Matrix realM = transformImageToMatrix(&realImage);
		FindBInA(realM, m);
	*/
	struct timeb start, end;
	Matrix M = MatReadAsc("matrice.mx");
	// Image IMMM = ImRead("Stereo/im0.pgm");

	ftime(&start);
	int C = WeylNorm(M);
	ftime(&end);
	int time_taken = (int)(1000.0 * (end.time - start.time) + (end.millitm - start.millitm));
	if (DEBUG)
		printf("Weyl:\n\tNorme attendue=193917\n\tNorme obtenue=%d\n\tTime taken: %dms\n\n", C, time_taken);

	ftime(&start);
	int D = WeylNormOpti(M);
	ftime(&end);
	time_taken = (int)(1000.0 * (end.time - start.time) + (end.millitm - start.millitm));
	if (DEBUG)
		printf("Weyl Opti:\n\tNorme attendue=193917\n\tNorme obtenue=%d\n\tTime taken: %dms\n", D, time_taken);

	Image IM0 = ImRead("Stereo/imB1.pgm");
	Matrix im0 = transformImageToMatrix(&IM0);
	Image IM1 = ImRead("Stereo/imB0.pgm");
	Matrix im1 = transformImageToMatrix(&IM1);

	ftime(&start);
	Matrix stereoResult = Stereo(im0, im1, 5, atoi(argv[1]), atoi(argv[2]));
	ftime(&end);
	time_taken = (int)(1000.0 * (end.time - start.time) + (end.millitm - start.millitm));
	printf("Stereo success\n");
	if (DEBUG)
		printf("Stereo:\n\tMatrix [%d %d]\n\tTime taken: %dms\n", D, ImNbRow(im0),ImNbCol(im0),time_taken);
	MatWriteAsc(stereoResult, "Stereo/dispMap.mx");

	// Matrix neigh1 = GetNeighborhood(m, 5, 5, 1);
	// MatWriteAsc(neigh1, "neigh1.mx");
	// Matrix neigh3 = GetNeighborhood(m, 5, 5, 3);
	// MatWriteAsc(neigh3, "neigh3.mx");
}