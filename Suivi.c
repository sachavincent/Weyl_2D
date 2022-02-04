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

void FindBInANeighborhood(Matrix A, Matrix B, int *posX, int *posY, int minFound)
{
    int prevX = *posX;
    int prevY = *posY;

    int heightA = MatNbRow(A);
    int widthA = MatNbCol(A);

    int heightB = MatNbRow(B);
    int widthB = MatNbCol(B);
    if (heightB > heightA || widthB > widthA)
    {
        printf("Error in findBInA, image B is larger than A\n");
        return;
    }
    // int **aValues = MatGetInt(A);
    int over = 1;
    int halfNeigh = 1;
    printf("Looking at (%d, %d)\n", *posX, *posY);

    for (int x = *posX - halfNeigh; x <= *posX + halfNeigh; x++)
    {
        for (int y = *posY - halfNeigh; y <= *posY + halfNeigh; y++)
        {
            Matrix differenceMatrix = GetSlicedDifferenceMatrix(A, x - widthB / 2 + 1, y - heightB / 2 + 1, B);

            int newNorm = WeylNorm(differenceMatrix);
            if (newNorm < minFound)
            {
                minFound = newNorm;
                if (x != *posX || y != *posY)
                {
                    over = 0;
                    *posX = x;
                    *posY = y;
                }
                printf("Better norm (= %d) found at (%d, %d)\n", minFound, x, y);
            }
            MatFree(&differenceMatrix);
        }
    }

    if (!over)
    {
        printf("For pixel (%d, %d), found best at (%d, %d) min=%d\n", prevX, prevY, *posX, *posY, minFound);
        FindBInANeighborhood(A, B, posX, posY, minFound);
    }
    else
        printf("Stopped at (%d, %d)\n", *posX, *posY);
}

void Suivi(Matrix imagette, int nbFrames)
{
    //int x = 302 + 132 / 2;
    //int y = 214 + 430 / 2;
    int x = 106;
    int y = 295;
    int heightB = MatNbRow(imagette);
    int widthB = MatNbCol(imagette);
    for (int i = 0; i < nbFrames; i++)
    {
        // char *fileName = malloc(19 + 5 + 4);
        // sprintf(fileName, "%s%d%s", "Suivi/Bike/img/frame_", (i + 1), ".pgm");
        char *fileName = malloc(19);
        sprintf(fileName, "%s%d%s", "Suivi/frame_", (i + 1), ".pgm");
        Matrix im_i = ReadMatrixFromPGM(fileName);
        if (im_i == NULL)
        {
            printf("Stopped at frame number %d because '%s' was not found!\n", i, fileName);
            break;
        }

        FindBInANeighborhood(im_i, imagette, &x, &y, INT_MAX);
        int heightIm = MatNbRow(im_i);
        int widthIm = MatNbCol(im_i);

        int **imValues = MatGetInt(im_i);

#pragma omp for
        for (int Vx = max(0, x - (widthB / 2)); Vx <= min(widthIm - 1, x + (widthB / 2)); Vx++)
        {
            for (int Vy = max(0, y - (heightB / 2)); Vy <= min(heightIm - 1, y + (heightB / 2)); Vy++)
            {
                if (Vx == x - (widthB / 2) || Vx == x + (widthB / 2) || Vy == y - (heightB / 2) || Vy == y + (heightB / 2))
                {
                    imValues[Vy][Vx] = 255;
                }
            }
        }

        Image image = TransformMatrixToImage(&im_i);
        // char *newFileName = malloc(25 + 5);
        // sprintf(newFileName, "%s%d%s", "Suivi/Bike/suivi_frame_", (i + 1), ".pgm");
        char *newFileName = malloc(25);
        sprintf(newFileName, "%s%d%s", "Suivi/suivi_frame_", (i + 1), ".pgm");
        printf("Successfully tracked frame %d at (%d, %d) saved at %s\n", i + 1, x, y, newFileName);
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

    // char *nomImagettePath = malloc(strlen(nomImagette) + 5 + 6);
    // sprintf(nomImagettePath, "%s%s", "Suivi/Bike/", nomImagette);
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