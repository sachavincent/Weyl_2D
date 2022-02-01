#include <assert.h>
#include <sys\timeb.h>
#include <errno.h>
#include <limits.h>

#include "utils.c"

char *directory = "Suivi";

int DEFAULT_THRESHOLD = 0;
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

void FindBInA(Matrix A, Matrix B, int *posX, int *posY)
{
    int heightA = MatNbRow(A);
    int widthA = MatNbCol(A);

    int heightB = MatNbRow(B);
    int widthB = MatNbCol(B);
    if (heightB > heightA || widthB > widthA)
    {
        printf("Error in findBInA, image B is larger than A\n");
        return;
    }
    int **aValues = MatGetInt(A);
    int normMin = INT_MAX;
    int over = 1;
    printf("B3\n");

    for (int x = *posX - 1; x <= *posX + 1; ++x)
    {
        for (int y = *posY - 1; y <= *posY + 1; ++y)
        {
            if (x == *posX && y == *posY)
                continue;
            Matrix differenceMatrix = GetSlicedDifferenceMatrix(A, x, y, B);
            int newNorm = WeylNorm(differenceMatrix);
            if (newNorm < normMin)
            {
                over = 0;
                *posX = x;
                *posY = y;
                normMin = newNorm;
            }
            MatFree(&differenceMatrix);
        }
    }

    printf("B6\n");

    if (over || *posX <= widthB / 2 || *posY <= heightB / 2)
    {
        for (int x = *posX - (widthB / 2); x <= *posX + (widthB / 2); ++x)
        {
            for (int y = *posY - (heightB / 2); y <= *posY + (heightB / 2); ++y)
            {
                if (x == *posX - (widthB / 2) || x == *posX + (widthB / 2) || y == *posY - (heightB / 2) || y == *posY + (heightB / 2))
                    aValues[y][x] = 255;
            }
        }
        return;
    }
    printf("B7\n");
    FindBInA(A, B, posX, posY);
    printf("B8\n");
}

void Suivi(Matrix imagette, int nbFrames)
{
    for (int i = 0; i < nbFrames; i++)
    {
        char *fileName = malloc(19);
        sprintf(fileName, "%s%d%s", "Suivi/frame_", (i + 1), ".pgm");
        Matrix im_i = ReadMatrixFromPGM(fileName);
        if (im_i == NULL)
        {
            printf("Stopped at frame number %d because '%s' was not found!", i, fileName);
            break;
        }

        printf("B4\n");
        int x = 105;
        int y = 295;
        FindBInA(im_i, imagette, &x, &y);
        printf("B5\n");

        Image image = TransformMatrixToImage(&im_i);
        char *newFileName = malloc(25);
        sprintf(newFileName, "%s%d%s", "Suivi/suivi_frame_", (i + 1), ".pgm");
        printf("Successfully tracked frame %d at %s", i + 1, newFileName);
        ImWriteAsc(image, newFileName);

        MatFree(&im_i);
        ImFree(&image);
        
        free(fileName);
        free(newFileName);
    }
}

int CheckImagette(char *imgArg)
{
    size_t len = strlen(imgArg);
    if (len < 5)
        return 1;

    char *ext = strrchr(imgArg, '.');
    if (!ext)
    {
        return 1;
    }

    return strcmp(".pgm", &imgArg[len - 4]);
}

int CheckNbFrames(char *nbFramesArg)
{
    char *end;
    long nbFramesLong = strtol(nbFramesArg, &end, 10);
    if (end == nbFramesArg || *end != '\0' || errno == ERANGE)
    {
        printf("Usage: ./Suivi [--opti] <imagette.pgm> <N frames>\n");
        exit(1);
    }

    int nbFrames = (int)nbFramesLong;

    if (nbFrames <= 0)
    {
        printf("N frames should be above 1\n");
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

    if ((opti && argc != 4) || (!opti && argc != 3))
    {
        printf("Usage: ./Suivi [--opti] <imagette.pgm> <N frames>\n");
        return 1;
    }

    char *nomImagette = opti ? argv[2] : argv[1];

    if (CheckImagette(nomImagette))
    {
        printf("Incorrect imagette name: %s!\n", nomImagette);
        return 3;
    }

    int nFrames = CheckNbFrames(opti ? argv[3] : argv[2]);

    char *nomImagettePath = malloc(strlen(nomImagette) + 6);
    sprintf(nomImagettePath, "%s%s", "Suivi/", nomImagette);

    printf("Imagette name: %s!\n", nomImagettePath);
    Matrix boundingBoxPGM = ReadMatrixFromPGM(nomImagettePath);
    if (boundingBoxPGM == NULL)
    {
        printf("Imagette not found: %s!\n", nomImagettePath);
        return 4;
    }

    Suivi(boundingBoxPGM, nFrames);

    MatFree(&boundingBoxPGM);
    free(nomImagettePath);

    return 0;
}