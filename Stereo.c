#include "utils.h"

char *directory = "Stereo";

int DEFAULT_NEIGH = 3;

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
 * @brief Stereovision
 *
 * @param Im1 Image on the left
 * @param Im2 Image on the right
 * @param neigh neighborhood
 * @param opti if 1 then use Weyl Opti
 * @param bidirectional if 1 then invert Stereo
 *
 * @return Matrix disparity Map
 */
Matrix Stereo(Matrix Im1, Matrix Im2, int neigh, int opti, int bidirectional)
{
    assert(neigh % 2 == 1);

    int width1 = MatNbCol(Im1);
    int height1 = MatNbRow(Im1);
    int width2 = MatNbCol(Im2);
    int height2 = MatNbRow(Im2);

    assert(width1 > 0 && height1 > 0);
    assert(width1 == width2 && height1 == height2);

    int halfNeigh = neigh / 2;

    Matrix disparityMap = MatAlloc(Int, height1 - halfNeigh * 2, width1 - halfNeigh * 2);
    int **mapValues = MatGetInt(disparityMap);

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
            int distance = 0;
            Matrix imagette = GetNeighborhood(Im1, x, y, neigh);

            int startSearch = bidirectional ? y : halfNeigh;
            int endSearch = bidirectional ? width1 - halfNeigh : y + 1;
            // Recherche dans Im2
            for (int y2 = startSearch; y2 < endSearch; y2++)
            {
                Matrix diff;
                if (bidirectional)
                    diff = DifferenceMatrices(listImagettes[y2 - halfNeigh], imagette);
                else
                    diff = DifferenceMatrices(imagette, listImagettes[y2 - halfNeigh]);
                int norm;
                if (opti)
                    norm = WeylNormOpti(diff);
                else
                    norm = WeylNorm(diff);

                if (norm <= bestNorm)
                {
                    bestNorm = norm;
                    distance = abs(y2 - y);
                }

                MatFree(&diff);
            }
            MatFree(&imagette);
            mapValues[x - halfNeigh][y - halfNeigh] = distance;
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

    MatFree(&toGetSize);

    return disparityMap;
}

int CheckNeigh(char *neighArg)
{
    char *end;
    long neighLong = strtol(neighArg, &end, 10);
    if (end == neighArg || *end != '\0' || errno == ERANGE)
    {
        printf("Usage: ./Stereo [--opti] [--bidirectional] <Image Left> <Image Right> [Neighborhood (default=%d)]\n", DEFAULT_NEIGH);
        exit(1);
    }

    int neigh = (int)neighLong;

    if (neigh % 2 == 0)
    {
        printf("Neighborhood should be an odd number\n");
        exit(1);
    }

    if (neigh <= 0 || neigh > 11)
    {
        printf("Neighborhood should be between 1 and 11\n");
        exit(1);
    }

    return neigh;
}

int main(int argc, char *argv[])
{
    int c;
    const char *short_opt = "b";
    struct option long_opt[] = {
        {"opti", no_argument, NULL, 'o'}};

    int opti = 0;
    int bidirectional = 0;
    while ((c = getopt_long(argc, argv, short_opt, long_opt, NULL)) != -1)
    {
        switch (c)
        {
        case -1: /* no more arguments */
        case 0:  /* long options toggles */
            break;

        case 'o':
            opti = 1;
            break;
        case 'b':
            bidirectional = 1;
            break;

        default:
            fprintf(stderr, "%s: invalid option -- %c\n", argv[0], c);
            fprintf(stderr, "Try `%s --help' for more information.\n", argv[0]);
            return 2;
        };
    };
    int idxNeigh = 3 + opti + bidirectional;

    if (argc != idxNeigh + 1 && argc != idxNeigh)
    {
        printf("Usage: ./Stereo [--opti] [--bidirectional] <Image Left> <Image Right> [Neighborhood (default=%d)]\n", DEFAULT_NEIGH);
        return 1;
    }
    int neigh = DEFAULT_NEIGH;
    if (argc >= idxNeigh + 1)
        neigh = CheckNeigh(argv[idxNeigh]);

    struct timeb start, end, startGlobal, endGlobal;

    ftime(&startGlobal);
    char *imageLeft = argv[1 + opti + bidirectional];
    char *imageRight = argv[2 + opti + bidirectional];

    char *pathLeft = malloc(strlen(imageLeft) + strlen(directory) + 2);
    sprintf(pathLeft, "%s%s%s", directory, "/", imageLeft);
    char *pathRight = malloc(strlen(imageLeft) + strlen(directory) + 2);
    sprintf(pathRight, "%s%s%s", directory, "/", imageRight);

    if (DEBUG)
    {
        printf("pathLeft: %s\n", pathLeft);
        printf("pathRight: %s\n", pathRight);
    }

    Matrix matImLeft = ReadMatrixFromPGM(pathLeft);
    Matrix matImRight = ReadMatrixFromPGM(pathRight);

    ftime(&start);
    Matrix stereoResult = Stereo(matImLeft, matImRight, neigh, opti, 0);
    ftime(&end);

    int time_taken = (int)(1000.0 * (end.time - start.time) + (end.millitm - start.millitm));
    int time_takenB = 0;

    char neighPath[3];
    sprintf(neighPath, "%d", neigh);
    char *optiPath = opti ? "_OPTI" : "";

    char *finalPath;
    if (bidirectional)
    {
        char *bidirBPath = bidirectional ? "_B" : "";
        char *pathB = (char *)malloc(strlen(directory) + 11 + sizeof(neighPath) + sizeof(optiPath) + sizeof(bidirBPath) + 3 + 1);

        sprintf(pathB, "%s%s%s%s%s%s", directory, "/disparity_", neighPath, optiPath, bidirBPath, ".mx");

        ftime(&start);
        Matrix stereoResult2 = Stereo(matImRight, matImLeft, neigh, opti, 1);
        ftime(&end);
        time_takenB = (int)(1000.0 * (end.time - start.time) + (end.millitm - start.millitm));

        if (DEBUG)
            MatWriteAsc(stereoResult2, pathB);

        int height = MatNbRow(stereoResult);
        int width = MatNbCol(stereoResult);
        int **values = MatGetInt(stereoResult);
        int **values2 = MatGetInt(stereoResult2);

        Matrix finalMatrix = MatAlloc(Int, height, width);
        int **finalMatrixValues = MatGetInt(finalMatrix);
        for (int x = 0; x < height; x++)
        {
            for (int y = 0; y < width; y++)
            {
                finalMatrixValues[x][y] = values[x][y] == values2[x][y] ? values[x][y] : 0;
            }
        }

        finalPath = (char *)malloc(strlen(directory) + 11 + sizeof(neighPath) + sizeof(optiPath) + 6 + 3 + 1);

        sprintf(finalPath, "%s%s%s%s%s%s", directory, "/disparity_", neighPath, optiPath, "_BIDIR", ".mx");
        MatWriteAsc(finalMatrix, finalPath);
        MatFree(&finalMatrix);
        MatFree(&stereoResult2);

        free(pathB);
    }

    char *bidirAPath = bidirectional ? "_A" : "";

    char *pathA = (char *)malloc(strlen(directory) + 11 + sizeof(neighPath) + sizeof(optiPath) + sizeof(bidirAPath) + 3 + 1);

    sprintf(pathA, "%s%s%s%s%s%s", directory, "/disparity_", neighPath, optiPath, bidirAPath, ".mx");

    if (!bidirectional || DEBUG)
        MatWriteAsc(stereoResult, pathA);

    if (!bidirectional)
    {
        finalPath = (char *)malloc(strlen(directory) + 11 + sizeof(neighPath) + sizeof(optiPath) + 3 + 1);

        sprintf(finalPath, "%s%s%s%s%s", directory, "/disparity_", neighPath, optiPath, ".mx");
    }

    ftime(&endGlobal);
    int time_taken_global = (int)(1000.0 * (endGlobal.time - startGlobal.time) + (endGlobal.millitm - startGlobal.millitm));
    if (bidirectional)
    {
        printf("Stereo%ssuccessful:\n\tNeighborhood=%dx%d\n\tMatrix: [%d %d]\n\tTime taken: %d ms\n\t\tStereo L->R: %d ms\n\t\tStereo R->L: %d ms\n\tSaved at: %s\n", opti ? " (optimized) " : " ", neigh, neigh, MatNbRow(stereoResult), MatNbCol(stereoResult), time_taken_global, time_taken, time_takenB, finalPath);
    }
    else
        printf("Stereo%ssuccessful:\n\tNeighborhood=%dx%d\n\tMatrix: [%d %d]\n\tTime taken: %d ms\n\t\tStereo L->R: %d ms\n\tSaved at: %s\n", opti ? " (optimized) " : " ", neigh, neigh, MatNbRow(stereoResult), MatNbCol(stereoResult), time_taken_global, time_taken, finalPath);

    free(pathLeft);
    free(pathRight);
    free(pathA);
    free(finalPath);

    MatFree(&stereoResult);
    MatFree(&matImLeft);
    MatFree(&matImRight);
}