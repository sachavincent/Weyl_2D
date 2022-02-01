Matrix Stereo2(Matrix Im1, Matrix Im2, int neigh)
{
	int width1 = MatNbCol(Im1);
	int height1 = MatNbRow(Im1);
	int width2 = MatNbCol(Im2);
	int height2 = MatNbRow(Im2);

	assert(width1 > 0 && height1 > 0);
	assert(width1 == width2 && height1 == height2);

	int halfNeigh = neigh / 2;

	Matrix disparityMap = MatAlloc(Int, height1, width1);
	int **mapValues = MatGetInt(disparityMap);

	int maxNorm = INT_MIN;
	int maxValue = INT_MIN;

	Matrix toGetSize = MatAlloc(Int, neigh, neigh);

	// Parcours de Im1
#pragma omp for
	for (int x = halfNeigh; x < height1 - halfNeigh; x++)
	{
		//	printf("On est a x = %d ...\n", x);
		Matrix *listImagettes1 = (Matrix *)malloc(sizeof(toGetSize) * (width1 - halfNeigh * 2));
        Matrix *listImagettes2 = (Matrix *)malloc(sizeof(toGetSize) * (width1 - halfNeigh * 2));

		for (int y = halfNeigh; y < width1 - halfNeigh; y++)
		{
            listImagettes1[y - halfNeigh] = GetNeighborhood(Im2, x, y, neigh);
			listImagettes2[y - halfNeigh] = GetNeighborhood(Im2, x, y, neigh);
		}
		for (int y = halfNeigh; y < width1 - halfNeigh; y++)
		{
			int bestY;
			int bestNorm = INT_MAX;
			for (int y2 = y; y2 < width1 - halfNeigh; y2++)
			{
				int norm = norms[y2 - halfNeigh];
				if (norm < bestNorm )//&& bestNorm < 82609)
				{
					bestNorm = norm;
					bestY = y2;
				}
				if (maxNorm < norm)
				{
					maxNorm = norm;
				}
			}
			/*if (y == 340 && x == 475)
			{
				printf("Best norm=%d\n", bestNorm);
			}*/
			mapValues[x][y] = (int)(255.0 * (double)abs(bestY - y) / (double)width1);
			if (maxValue < mapValues[x][y])
				maxValue = mapValues[x][y];
		}
	}
	#pragma omp barrier
	return disparityMap;
}