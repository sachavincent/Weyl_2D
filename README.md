Implementation of 2D Weyl Norm

Two use cases :
- Stereovision : creates a disparity map using this norm
  This dataset can be used: https://vision.middlebury.edu/stereo/data/
- Tracking : tries to track an object frame by frame
  This dataset can be used: https://amoudgl.github.io/tlp/
  
Matlab/Octave scripts are included to help :
- "saveToPGM.m" : convert any image format to ".pgm" (Grayscale)
- "convertDatasetToPGM.m" : convert an entire Dataset to ".pgm"
- "convertVideoToFrames.m" : convert video to independant ".pgm" images
- "loadMXFile.m" : read ".mx" file (LiMaCe) as matrix
- "convertTrackingToVideo.m" : convert tracked frames to an ".avi" video format
