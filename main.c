#include <stdio.h>
#include <string.h>
#include "limace.h"

#define min(x, M) ((x) >= M ? M : x)

double WeylDistance(Matrix mat)
{
	return 0;
}
Matrix GetIntegralImage(Matrix m)
{
	// nbLignes = nx
	// nbColonnes = ny
	int height = MatNbRow(m);
	int width = MatNbCol(m);
	int **i = MatGetInt(m);

	Matrix mm = MatAlloc(Int, width, height);
	int **ii = MatGetInt(mm);

	//printf("%d\n", i[0][0]);
	
	ii[0][0] = i[0][0];
	for (int x = 1; x < height; x++)
	{
		ii[x][0] += ii[x - 1][0] + i[x][0];
	}

	for (int y = 1; y < width; y++)
	{
		int s = i[0][y]; // scalar accumulator
		ii[0][y] = ii[0][y - 1] + s;
		for (int x = 1; x < height; x++)
		{
			s += i[x][y];
			ii[x][y] = ii[x][y - 1] + s;
		}
	}
	printf("Success\n");
	/*
	int value = 0;
	for (int y = 0; y < width; ++y)
	{
		for (int x = 0; x < height; ++x) // for(unsigned int x = 0; x < WIDTH; ++x) {
		{

			value = i[x][y];
			if (y > 0)
				value += ii[x][y - 1];
			if (x > 0)
				value += ii[x - 1][y];
			if (x > 0 && y > 0)
			{
				value -= ii[x - 1][y - 1];
				//printf("%d = %d + %d + %d - %d \n",value,i[x][y],i[x][y-1],i[x-1][y],i[x-1][y-1]);
			}
			ii[x][y] = value;
		}
	}*/
	return mm;
}
/*
Image GetIntegralImage(Image m)
{
	// nbLignes = nx
	// nbColonnes = ny
	int height = ImNbRow(m);
	int width = ImNbCol(m);
	int **i = ImGetI(m);

	Image mm = ImAlloc(GrayLevel, width, height);
	unsigned char **ii = ImGetI(mm);

	//printf("%d\n", i[0][0]);
	int value = 0;
	for (int y = 0; y < width; ++y)
	{
		for (int x = 0; x < height; ++x) // for(unsigned int x = 0; x < WIDTH; ++x) {
		{
			
       		value = i[x][y];
			if (y > 0)
				value += ii[x][y - 1];
			if (x > 0)
				value += ii[x - 1][y];
			if (x > 0 && y > 0){
				value -= ii[x - 1][y - 1];			
				//printf("%d = %d + %d + %d - %d \n",value,i[x][y],i[x][y-1],i[x-1][y],i[x-1][y-1]);
			}
			ii[x][y] = value;
		}
	}
	return mm;
}
*/

/*
int getIntegralImage(Matrix m, int x0, int y0)
{
	int h = MatNbRow(m);
	int w = MatNbCol(m);
	
	//Matrix M_TL = getIntegralPart(TL);
	//Matrix M_TR = getIntegralPart(TR);
	//Matrix M_BL = getIntegralPart(BL);
	//Matrix M_BR = getIntegralPart(BR);

	int **ii = MatGetInt(m);
	// Reconstruction de l'image integrale :
	// Computing the sum within a rectangle
	int A = 0, B = 0, C = 0, D = 0;
	if (x0 > 0 && y0 > 0)
		A = ii[x0 - 1][y0 - 1];
	if (y0 > 0)
		B = ii[x0 + w - 1][y0 - 1];
	if (x0 > 0)
		C = ii[x0 - 1][y0 + h - 1];

	D = ii[x0 + w - 1][y0 + h - 1];

	int sum = A - B - C + D;

	return sum;
}
Image getIntegralImage(Matrix m)
{
	return getIntegralImage(m, 0, 0);
}
*/

Matrix getIntegralPart(Matrix m)
{
	// nbLignes = nx
	// nbColonnes = ny
	int nx = MatNbRow(m);
	int ny = MatNbCol(m);
	Matrix mm = MatAlloc(Int, nx, ny);

	int **i = MatGetInt(m);
	int **ii = MatGetInt(mm);

	ii[0][0] = i[0][0];
	//printf("%d\n", i[0][0]);
	for (int x = 1; x < nx; x++)
	{
		ii[x][0] += ii[x - 1][0] + i[x][0];
	}

	for (int y = 1; y < ny; y++)
	{
		int s = i[0][y]; // scalar accumulator
		ii[0][y] = ii[0][y - 1] + s;
		for (int x = 1; x < nx; x++)
		{
			s += i[x][y];
			ii[x][y] = ii[x][y - 1] + s;
		}
	}

	return mm;
}

int main(int argc, char *argv[])
{
	if (argc != 2)
	{
		printf("Error parameter");
		return 1;
	}

	//Image m = ImRead(argv[1]);
	Matrix m = MatReadAsc(argv[1]);
	if (m == NULL)
	{
		printf("Probleme lors de la lecture de la matrice");
		return 2;
	}
	//Image IG = GetIntegralImage(m);
	Matrix M_IG = GetIntegralImage(m);
	//ImWriteAsc(IG, "test.mx");
	MatWriteAsc(M_IG, "test.mx");

	//	Matrix integralImage = getIntegralImage(m);
	//	double dist = WeylDistance(integralImage);
	// Norme attendue = 193917
}