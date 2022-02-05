function convertTrackingToVideo(directory, nbFrames, name);
pkg load video;

listFiles = dir(directory);
if nargin < 2
  listFiles = listFiles(3:end);
  nbFrames = rows(listFiles);
else
  listFiles = listFiles(3:nbFrames+3);
endif

if nargin < 3
  name = 'tracking_result.avi';
endif
nextPrint = 10;
resVid = VideoWriter(strcat(directory, '/', name));
frameName = strcat(directory, "/tracking_frame_1.pgm");

fig = figure('Toolbar','none','Menubar','none');
ax = gca;
set(ax, 'Position', [0 0 1 1]);
set(fig, 'Visible', 'off');

if !isfile(frameName)
  disp("Couldn't find any tracking files.");
  return;
endif

disp(["Converting ", num2str(nbFrames), " frames..."]);
for i = 1:nbFrames
  frameName = strcat(directory, "/tracking_frame_", int2str(i), ".pgm");
  if !isfile(frameName)
    disp(["Stopped convertion at frame", int2str(i-1),"!"]);
    break;
  endif
  
  im = imread(frameName);
  
  set(fig, 'Position', [0 0 columns(im) rows(im)]);
  imshow(im);
  drawnow
  writeVideo (resVid, getframe(gcf));
  
  if i*100/nbFrames >= nextPrint
    disp(["Conversion is ", num2str(nextPrint), "% done..."]);
    nextPrint += 10;
  endif
endfor
close(resVid);

disp(["Conversion successful!"]);

figure;
for i=1:71
  imshow(['Suivi\University\frames\frame_', int2str(i), '.jpg']);
  drawnow
endfor