pkg load video;
vid = VideoReader ('Suivi/video.webm');
numFrames = vid.NumberOfFrames;
n = numFrames;
im = [];
i = 1;
while (!isempty (img = readFrame(vid)))
    imwrite(img,['Suivi/frame_' int2str(i), '.jpg']);
    i++;
    break;
endwhile
close(vid);


vid = VideoReader ('Suivi/video.webm');
resVid = VideoWriter('Suivi/videoRes.avi');
%resVid.FrameRate = vid.FrameRate;

numFrames = vid.NumberOfFrames;
im = [];
numFrames = 71;
for i = 1:numFrames
  figure;
  imshow(imread(["Suivi/suivi_frame_", int2str(i), ".pgm"]));
  drawnow
  ax = gca;
  set(ax,"Units",'pixels');
  pos = get(ax, "Position");
  ti = get(ax, "TightInset");
  rect = [ceil(pos(1)), ceil(pos(2)+ti(2)/2)+1, pos(3)-1, floor(pos(4)-1-ti(2))];
  F = getframe(gca,rect);
  set(ax, "Units",'normalized');
  writeVideo (resVid, F.cdata);
  close;
endfor
close(resVid);
close(vid);

##############################
frame_1 = imread("Suivi/suivi_frame_1.pgm");
figure;
imshow(frame_1);
figure;
frame_2 = imread("Suivi/suivi_frame_2.pgm");
imshow(frame_2);
figure;
frame_5 = imread("Suivi/suivi_frame_10.pgm");
imshow(frame_5);

figure;
imshow(outputimg);
drawnow
imshow(F.cdata)
imwrite(F.cdata, "t.jpg")