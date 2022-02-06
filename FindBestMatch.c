#include "utils.h"

char *directory = "BestMatch";

/**
 * @brief Find best match of B in A
 *
 * @param A full matrix
 * @param B smaller matrix
 * @param opti if 1 then use Weyl Opti
 * @param disparityMatrix if 1 then save disparity map
 *
 * @return Matrix best match for B in A
 */
Matrix FindBestMatch(Matrix A, Matrix B, int opti, int saveDispMap)
{
    int widthA = MatNbCol(A);
    int heightA = MatNbRow(A);
    int widthB = MatNbCol(B);
    int heightB = MatNbRow(B);

    assert(widthA > 0 && heightA > 0);
    assert(widthB > 0 && heightB > 0);
    assert(widthA >= widthB || heightA >= heightB);

    Matrix bestMatch = MatAlloc(Int, heightB, widthB);

    int normMin = INT_MAX;

    Matrix dispMatrix = NULL;
    int **dispMap;
    if (saveDispMap)
    {
        dispMatrix = MatAlloc(Int, heightA - heightB, widthA - widthB);
        dispMap = MatGetInt(dispMatrix);
    }

    if (DEBUG)
    {
        printf("Matrix A of size %dx%d\n", widthA, heightA);
        printf("Matrix B of size %dx%d\n", widthB, heightB);

        printf("maxX = %d, maxY = %d\n", heightA - heightB, widthA - widthB);
    }

    for (int x = 0; x < heightA - heightB; x++)
    {
        for (int y = 0; y < widthA - widthB; y++)
        {
            Matrix differenceMatrix = GetSlicedDifferenceMatrix(A, y, x, B);
            int newNorm;
            if (opti)
                newNorm = WeylNormOpti(differenceMatrix);
            else
                newNorm = WeylNorm(differenceMatrix);

            if (saveDispMap)
                dispMap[x][y] = newNorm;

            if (newNorm < normMin)
            {
                normMin = newNorm;
                bestMatch = SliceMatrix(A, y, x, widthB, heightB);
            }
            MatFree(&differenceMatrix);
        }
    }
#pragma omp barrier
    if (DEBUG)
    {
        printf("FindBestMatch done => minNorm=%d\n", normMin);
    }

    if (saveDispMap)
    {
        char *dispFileName = "disparity_map.mx";
        char *pathDisp = malloc(strlen(directory) + strlen(dispFileName) + 2);
        sprintf(pathDisp, "%s%s%s", directory, "/", dispFileName);

        MatWriteAsc(dispMatrix, pathDisp);
        MatFree(&dispMatrix);
        free(pathDisp);
    }

    return bestMatch;
}

int main(int argc, char *argv[])
{
    int c;
    const char *short_opt = "";
    struct option long_opt[] = {
        {"opti", no_argument, NULL, 'o'}};

    int opti = 0;
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

        default:
            fprintf(stderr, "%s: invalid option -- %c\n", argv[0], c);
            fprintf(stderr, "Try `%s --help' for more information.\n", argv[0]);
            return 2;
        };
    };
    if (argc < opti + 3 || argc > 4 + opti)
    {
        printf("Usage: ./FindBestMatch [--opti] <Image> <Pattern> [SaveDispMap (default=1)]\n");
        return 1;
    }

    char *nameImA = argv[1 + opti];
    char *nameImB = argv[2 + opti];
    int saveDispMap = 1;
    if (argc > 3 + opti)
        saveDispMap = atoi(argv[3 + opti]);

    struct timeb start, end;

    char *pathA = malloc(strlen(nameImA) + strlen(directory) + 2);
    sprintf(pathA, "%s%s%s", directory, "/", nameImA);
    char *pathB = malloc(strlen(nameImB) + strlen(directory) + 2);
    sprintf(pathB, "%s%s%s", directory, "/", nameImB);

    if (DEBUG)
    {
        printf("pathA: %s\n", pathA);
        printf("pathB: %s\n", pathB);
    }

    Matrix matrixA = ReadMatrixFromPGM(pathA);
    if (matrixA == NULL)
        return 1;
    Matrix matrixB = ReadMatrixFromPGM(pathB);
    if (matrixB == NULL)
    {
        free(pathA);
        free(pathB);
        MatFree(&matrixA);
        return 1;
    }

    ftime(&start);
    Matrix bestMatch = FindBestMatch(matrixA, matrixB, opti, saveDispMap);
    ftime(&end);
    int time_taken = (int)(1000.0 * (end.time - start.time) + (end.millitm - start.millitm));

    if (bestMatch == NULL)
    {
        printf("Error while trying to find best match...\n");
    }
    else
    {
        Image bestMatchImg = TransformMatrixToImage(&bestMatch);

        char *dispFileName = "best_match.pgm";
        char *pathBestMatch = malloc(strlen(directory) + strlen(dispFileName) + 2);
        sprintf(pathBestMatch, "%s%s%s", directory, "/", dispFileName);

        ImWriteAsc(bestMatchImg, pathBestMatch);

        if (saveDispMap)
        {
            char *dispFileName = "disparity_map.mx";
            char *pathDisp = malloc(strlen(directory) + strlen(dispFileName) + 2);
            sprintf(pathDisp, "%s%s%s", directory, "/", dispFileName);

            printf("Successfully found best match%s, disparity map saved at '%s'. Time taken=%d ms\n", opti ? " (optimized)" : "", pathDisp, time_taken);
            free(pathDisp);
        }
        else
        {
            printf("Successfully found best match%s. Time taken=%d ms\n", opti ? " (optimized)" : "", time_taken);
        }
        MatFree(&bestMatch);
        ImFree(&bestMatchImg);
        free(pathBestMatch);
    }

    free(pathA);
    free(pathB);

    MatFree(&matrixA);
    MatFree(&matrixB);
}