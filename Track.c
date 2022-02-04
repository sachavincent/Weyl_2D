#include "utils.h"
#include <sys/stat.h>
#include <sys/types.h>

char *directory = "Suivi";

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
 * @brief Looks for B in A at given position using neighborhood, takes min as result
 *
 * @param A whole matrix
 * @param B bounding box
 * @param posX top left x coordinate
 * @param posY top left y coordinate
 * @param minFound min norm found
 * @param neigh neighborhood size
 * @param opti if 1 then use Weyl Opti
 */
void FindBInANeighborhood(Matrix A, Matrix B, int *posX, int *posY, int *minFound, const int neigh, const int opti)
{
    int prevX = *posX;
    int prevY = *posY;

    int heightA = MatNbRow(A);
    int widthA = MatNbCol(A);

    int heightB = MatNbRow(B);
    int widthB = MatNbCol(B);
    assert(heightB <= heightA || widthB <= widthA);

    int done = 1;
    int halfNeigh = neigh / 2;
    if (DEBUG)
        printf("Looking at (%d, %d)\n", *posX, *posY);

    int minX = max(0, *posX - halfNeigh);
    int minY = max(0, *posY - halfNeigh);
    int maxX = min(widthA - 1 - widthB, *posX + halfNeigh);
    int maxY = min(heightA - 1 - heightB, *posY + halfNeigh);
    for (int x = minX; x <= maxX; x++)
    {
        for (int y = minY; y <= maxY; y++)
        {
            if (x == prevX && y == prevY && *minFound != INT_MAX)
                continue;
            Matrix differenceMatrix = GetSlicedDifferenceMatrix(A, x, y, B);

            int newNorm;
            if (opti)
                newNorm = WeylNormOpti(differenceMatrix);
            else
                newNorm = WeylNorm(differenceMatrix);

            if (newNorm < *minFound)
            {
                *minFound = newNorm;
                done = 0;
                *posX = x;
                *posY = y;

                if (DEBUG)
                    printf("Better norm (= %d) found at (%d, %d)\n", *minFound, x, y);
            }
            MatFree(&differenceMatrix);
        }
    }
    if (prevX == *posX && prevY == *posY) // Best is at current position, we can stop
    {
        done = 1;
    }

    if (!done)
    {
        if (DEBUG)
            printf("For pixel (%d, %d), found best at (%d, %d) min=%d\n", prevX, prevY, *posX, *posY, *minFound);
        FindBInANeighborhood(A, B, posX, posY, minFound, neigh, opti);
    }
    else if (DEBUG)
        printf("Stopped at (%d, %d)\n", *posX, *posY);
}

/**
 * @param dataset name of the dataset (needs to be located at Suivi/)
 * @param x top left position of the Bounding Box on the first frame
 * @param y top left position of the Bounding Box on the first frame
 * @param width width of the Bounding Box
 * @param height height of the Bounding Box
 * @param nbFrames if 0 then use all frames available
 * @param neigh (odd) neighborhood size in which to look for bounding box in the next frame
 * @param opti if 1 then use Weyl Opti
 */
int Track(char *dataset, int x, int y, const int width, const int height, int nbFrames, const int neigh, const int opti)
{
    int startX = x;
    int startY = y;
    assert(neigh % 2 == 1);
    assert(nbFrames >= 0);
    assert(x >= 0 && y >= 0);

    nbFrames = nbFrames == 0 ? INT_MAX : nbFrames;
    Matrix boundingBox = NULL;

    char *trackingDir = malloc(sizeof(directory) + 1 + sizeof(dataset) + 1 + 17 + (opti ? 5 : 0));
    sprintf(trackingDir, "%s%s%s%s%s", directory, "/", dataset, "/tracking", opti ? "_opti" : "");
    mkdir(trackingDir);

    for (int i = 0; i < nbFrames; i++)
    {
        char *fileName = malloc(sizeof(directory) + 1 + sizeof(dataset) + 1 + 14);
        sprintf(fileName, "%s%s%s%s%s%d%s", directory, "/", dataset, "/", "frame_", (i + 1), ".pgm");
        Matrix im_i = ReadMatrixFromPGM(fileName);
        if (im_i == NULL)
        {
            printf("Stopped at frame number %d because '%s' was not found!\n", i + 1, fileName);
            free(fileName);
            break;
        }

        // Create Bounding Box on first frame
        if (i == 0)
        {
            boundingBox = SliceMatrix(im_i, startX, startY, width, height);
            if (boundingBox == NULL)
            {
                return 1;
            }

            MatWriteAsc(boundingBox, "tb.mx");
        }
        int minFound = INT_MAX;
        // Find best match on all the other frames
        FindBInANeighborhood(im_i, boundingBox, &x, &y, &minFound, neigh, opti);
        int heightIm = MatNbRow(im_i);
        int widthIm = MatNbCol(im_i);

        int **imValues = MatGetInt(im_i);

        // Draw bounding box at coordinates
#pragma omp for
        for (int Vx = max(0, x); Vx <= min(widthIm - 1, x + width); Vx++)
        {
            for (int Vy = max(0, y); Vy <= min(heightIm - 1, y + height); Vy++)
            {
                // if (Vx == x || Vx == (x + width) || Vy == y || Vy == (y + height))
                if ((Vx >= x - 1 && Vx <= x + 1) || (Vx >= x + width - 1 && Vx <= x + 1 + width) ||
                    (Vy >= y - 1 && Vy <= y + 1) || (Vy >= y + height - 1 && Vy <= y + 1 + height))
                {
                    imValues[Vy][Vx] = 255;
                }
            }
        }

        Image image = TransformMatrixToImage(&im_i);
        char *newFileName = malloc(sizeof(trackingDir) + 1 + 24);

        sprintf(newFileName, "%s%s%d%s", trackingDir, "/tracking_frame_", (i + 1), ".pgm");
        printf("Successfully tracked%sbounding box at (%d, %d) on frame %d => saved at '%s'\n", opti ? " (optimized) " : " ", x, y, i + 1, newFileName);
        ImWriteAsc(image, newFileName);
        ImFree(&image);

        MatFree(&im_i);

        free(newFileName);
        free(fileName);
    }
    free(trackingDir);
    MatFree(&boundingBox);

    return 0;
}

int CheckNbFrames(char *nbFramesArg)
{
    char *end;
    long nbFramesLong = strtol(nbFramesArg, &end, 10);
    if (end == nbFramesArg || *end != '\0' || errno == ERANGE)
    {
        printf("Usage: ./Track [--opti] <Dataset name> <x (TopLeft)> <y (TopLeft)> <width> <height> (Bounding Box on 1st frame) [N frames (default=0=all)] [Neighborhood (default=%d)]\n", DEFAULT_NEIGH);
        exit(1);
    }

    int nbFrames = (int)nbFramesLong;

    if (nbFrames < 0)
    {
        printf("N frames should be >= 0\n");
        exit(1);
    }
    return nbFrames;
}

int main(int argc, char *argv[])
{
    int c;
    const char *short_opt = "";
    struct option long_opt[] = {{"opti", no_argument, NULL, 'o'}};

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

    if (argc < opti + 6 || argc > opti + 8)
    {
        printf("Usage: ./Track [--opti] <Dataset name> <x (TopLeft)> <y (TopLeft)> <width> <height> (Bounding Box on 1st frame) [N frames (default=0=all)] [Neighborhood (default=%d)]\n", DEFAULT_NEIGH);
        return 1;
    }

    char *datasetName = argv[1 + opti];
    int topX = atoi(argv[2 + opti]);
    int topY = atoi(argv[3 + opti]);
    int width = atoi(argv[4 + opti]);
    int height = atoi(argv[5 + opti]);
    int nbFrames = 0;
    int neigh = DEFAULT_NEIGH;

    if (argc > opti + 6)
    {
        nbFrames = CheckNbFrames(argv[6 + opti]);
    }

    if (argc > opti + 7)
    {
        neigh = atoi(argv[7 + opti]);
    }

    if (DEBUG)
    {
        printf("datasetName: %s\n", datasetName);
        printf("topX: %d\n", topX);
        printf("topY: %d\n", topY);
        printf("width: %d\n", width);
        printf("height: %d\n", height);
        printf("nbFrames: %d\n", nbFrames);
        printf("neigh: %d\n", neigh);
    }

    char str[3];
    if (nbFrames > 0)
        sprintf(str, "%d", nbFrames);
    else
        sprintf(str, "%s", "all");

    printf("Starting tracking of '%s' for %s frame%s...\n", datasetName, str, nbFrames == 0 || nbFrames > 1 ? "s" : "");

    struct timeb start, end;

    ftime(&start);
    int res = Track(datasetName, topX, topY, width, height, nbFrames, neigh, opti);
    ftime(&end);
    int time_taken = (int)(1000.0 * (end.time - start.time) + (end.millitm - start.millitm));

    if (res)
        printf("Error while tracking '%s'!\n", datasetName);
    else
        printf("Successfully tracked '%s' for %s frame%s in %f s!\n", datasetName, str, nbFrames == 0 || nbFrames > 1 ? "s" : "", (double)time_taken / 1000.0);

    return res;
}