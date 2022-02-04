function convertDatasetToPGM(directory);
listFiles = dir(directory);
for i=1:rows(listFiles)-2
  file = listFiles(i+2);
  fln = strcat(file.folder, '\', file.name);
  
  saveToPGM(strcat(file.folder, '\', file.name), strcat(file.folder, '\', "frame_", int2str(i))); 
endfor