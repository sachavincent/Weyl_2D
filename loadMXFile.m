function mat = loadMXFile(fileName);
  
fileMat = double(textread(fileName));
sizeX = fileMat(2);
sizeY = fileMat(3);
fileMat = fileMat(4:end)';
if(max(fileMat) <= 255 && min(fileMat) >= 0)
  fileMat = uint8(fileMat);
endif
mat = reshape(fileMat, [sizeY, sizeX]);
mat = transpose(mat);