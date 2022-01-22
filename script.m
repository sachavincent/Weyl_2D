pkg load image;
close all; clear all;
function mat = loadMXFile(fileName)
  fileMat = double(textread(fileName));
  sizeX = fileMat(2);
  sizeY = fileMat(3);
  fileMat = fileMat(4:end)';
  if(max(fileMat) <= 255 && min(fileMat) >= 0)
    fileMat = uint8(fileMat);
  endif
  mat = reshape(fileMat, [sizeY, sizeX]);
  mat = transpose(mat);
endfunction

A2 = loadMXFile("matrice.mx");
%A2 = double(imread("lena.jpg"));
R2 = loadMXFile("test.mx");
S = integralImage(transpose(A2));
colormap gray;
subplot(131); imagesc(A2); title("Image");
subplot(132); imagesc(R2); title("Image integrale");
subplot(133); imagesc(S); title("Ref Image integrale");

function J = II (I)
  ## FIXME: Can this part be more vectorized?
  s = size (I);
  s1 = s + [1,2];
  J = zeros (s1);
  J(2,2:s(2)+1) = I(1,:);
  s21 = s(2)+1;
  for y = 3:s1(1)
    y1 = y-1;
    J(y,1) = J(y1,2);
    J(y,2:s21) = J(y1,1:s(2)) + J(y1,3:s1(2)) - J(y-2,2:s21) + I(y1,:) + I(y-2,:); 
    J(y,end) = J(y1,s21);
  endfor
endfunction

figure;
colormap gray; imagesc(A2);
figure;
colormap gray; imagesc(II(A2));
figure;
colormap gray; imagesc(integralImage(A2, "rotated"));


% Disparity maps
full = imread("allCoins.pgm");
figure;
colormap gray; imagesc(full);
DispMat = loadMXFile("disparityMatrix.mx");
colormap gray; imagesc(DispMat);

% testing
test = loadMXFile("testFindBinA.mx");
colormap gray; imshow(test);