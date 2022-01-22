#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "limace.h"
#include <assert.h>

#define min(x, M) ((x) >= M ? M : x)
#define max(x, M) ((x) < M ? M : x)
int distance(int a, int b)
{
	return abs(b - a);
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

	int sizeNeigh = neigh / 2;
	int nbW = 0;

	printf("%d %d\n", m[xC][yC], sizeNeigh);
	Matrix res = MatAlloc(Int, neigh, neigh);
	int **mRes = MatGetInt(res);
	for (int x = xC - sizeNeigh; x <= xC + sizeNeigh; x++)
	{
		int nbH = 0;
		for (int y = yC - sizeNeigh; y <= yC + sizeNeigh; y++)
		{
			int pixel = m[x][y];

			mRes[nbW][nbH] = pixel;
			nbH++;
		}
		nbW++;
	}

	return res;
}

void Stereo(Matrix Im1, Matrix Im2, int neigh)
{
	int width1 = MatNbCol(Im1);
	int height1 = MatNbRow(Im1);
	int width2 = MatNbCol(Im2);
	int height2 = MatNbRow(Im2);

	assert(width1 == width2 && height1 == height2);
	int **im1 = MatGetInt(Im1);
	int **im2 = MatGetInt(Im2);

	int sizeNeigh = neigh / 2;
	for (int x = sizeNeigh; x < width1 - sizeNeigh; x++)
	{
		for (int y = sizeNeigh; y < height1 - sizeNeigh; y++)
		{
			int pixel = im1[x][y];
			// Recherche dans Im2
			for (int x2 = x; x2 < width1 - sizeNeigh; x2++)
			{
				int pixel2 = im2[x2][y];
			}
		}
	}
}

Matrix GetIntegralImage(Matrix m)
{
	int height = MatNbRow(m);
	int width = MatNbCol(m);
	int **i = MatGetInt(m);

	Matrix mm = MatAlloc(Int, height, width);
	int **ii = MatGetInt(mm);

	ii[0][0] = i[0][0];
	for (int x = 1; x < height; x++)
	{
		ii[x][0] = ii[x - 1][0] + i[x][0];
	}

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
	int pi1, pi2, pi3, pi4;
	pi1 = maxPi1 - minPi1;
	pi2 = maxPi2 - minPi2;
	pi3 = maxPi3 - minPi3;
	pi4 = maxPi4 - minPi4;
	return max(max(pi1, pi2), max(pi4, pi3));
}

Matrix GetSlicedMatrix(Matrix A, int x, int y, int h, int w)
{
	Matrix sliced = MatAlloc(Int, h, w);
	int **a = MatGetInt(A);
	int **m = MatGetInt(sliced);
	for (int i = x; i < x + h; i++)
	{
		for (int j = y; j < y + w; j++)
		{
			m[i - x][j - y] = a[i][j];
		}
	}
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

	int normB = WeylNorm(B);
	int distanceMin = INT_MAX;
	Matrix ABis = MatCopy(A);
	Matrix found = GetSlicedMatrix(ABis, 0, 0, heightB, widthB);
	//int **foundValues = MatGetInt(found);
	Matrix disparityMatrix = MatAlloc(Int, heightA - heightB, widthA - widthB);
	int **dm = MatGetInt(disparityMatrix);
	//Matrix slicedA = GetSlicedMatrix(A, 0, 0, heightB, widthB);
	// TODO: Ballek des bords

	for (int x = 0; x < heightA - heightB; x++)
	{
		for (int y = 0; y < widthA - widthB; y++)
		{
			Matrix slicedA = GetSlicedMatrix(A, x, y, heightB, widthB);
			int newDistance = distance(WeylNorm(slicedA), normB);
			dm[x][y] = newDistance;
			if (newDistance <= distanceMin)
			{
				distanceMin = newDistance;
				MatWriteAsc(slicedA, "testFindBInA.mx");
			}
		}
	}
	MatWriteAsc(disparityMatrix, "disparityMatrix.mx");
	return found;
}

Matrix transformImageToMatrix(Image *i)
{
	int height = ImNbRow(*i);
	int width = ImNbCol(*i);
	Matrix m = MatAlloc(Int, height, width);
	unsigned char **mat = ImGetI(*i);
	int **mi = MatGetInt(m);
	for (int x = 0; x < height; x++)
	{
		for (int y = 0; y < width; y++)
		{
			mi[x][y] = mat[x][y];
		}
	}
	MatWriteAsc(m, "testTransform.mx");

	printf("Converted Image Success\n");
	return m;
}

int main(int argc, char *argv[])
{
	if (argc != 2)
	{
		printf("Error parameter");
		return 1;
	}

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
	//Matrix neigh1 = GetNeighborhood(m, 5, 5, 1);
	//MatWriteAsc(neigh1, "neigh1.mx");
	//Matrix neigh3 = GetNeighborhood(m, 5, 5, 3);
	//MatWriteAsc(neigh3, "neigh3.mx");
}