function convertDatasetToPGM(directory, del);  
pkg load mapping;

listFiles = dir(directory);
sortedList = extractfield(listFiles, 'name')';
sortedList = sort_nat(sortedList)(3:end);
for i=1:rows(sortedList)
  fileName = cell2mat(sortedList(i));
  fln = strcat(directory, '/', fileName);
  saveToPGM(fln, strcat(directory, '/', "frame_", int2str(i))); 
  
  if del
    delete(fln);
  endif
endfor

disp(["Conversion successful!"]);