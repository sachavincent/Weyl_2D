function convertVideoToFrames(videoName, directory, nbFrames);
pkg load video;

if nargin < 2
  disp("directory needs to be specified");
  return;
endif  
if nargin < 3
  nbFrames = intmax;
endif

vid = VideoReader(videoName);
i = 1;
mkdir(directory);

while (!isempty (img = readFrame(vid)) && i <= nbFrames)
    imwrite(img,strcat(directory, '/frame_', int2str(i), '.jpg'));
    i++;
endwhile

close(vid);

disp(["Conversion successful!"]);