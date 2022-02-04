function saveToPGM(imgName, newName);
im = imread(imgName);
imwrite(im, strcat(newName, ".pgm"));