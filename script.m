pkg load image;
close all; clear all;

A2 = loadMXFile("matrice.mx");
%A2 = double(imread("lena.jpg"));
R2 = loadMXFile("test.mx");
S = integralImage(A2);
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

function l = lerp(a,b,t)
  l= a+t*(b-a);
end
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
figure;
colormap gray; imshow(imgette);

% Stereo
imdisp0 = imread("Stereo/disp0.pgm");
imdisp1 = imread("Stereo/disp1.pgm");
imdisp0 = imdisp0*(255 / max(imdisp0(:)));
imdisp1 = imdisp1*(255 / max(imdisp1(:)));
figure;
colormap gray; imshow(imdisp0);
figure;
colormap gray; imshow(imdisp1);

stereoDisp = loadMXFile("Stereo/disparity_0_7.mx");
stereoDisp = uint8(255*(stereoDisp - min(stereoDisp(:))) / (max(stereoDisp(:)) - min(stereoDisp(:))));
figure;colormap gray; 
subplot(131); 
imshow(imread("Stereo/tsukubagauche.pgm"));
subplot(132); 
imshow(imread("Stereo/tsukubadroite.pgm"));
subplot(133); imshow(stereoDisp);



D1 = loadMXFile("Stereo/disparity_11.mx");
figure;colormap gray;imagesc(D1);

D2 = loadMXFile("Stereo/bidir.mx");
figure;colormap gray;imagesc(D2);

D3 = loadMXFile("Stereo/disparity_11_BIDIR.mx");
D3 = double(D3);
D3 = uint8(255*((D3 - min(D3(:))) / (max(D3(:)) - min(D3(:)))));
figure;colormap gray;imshow(D3);