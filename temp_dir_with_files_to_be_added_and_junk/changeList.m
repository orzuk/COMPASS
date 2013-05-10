function changeList(dirName,fileName)
% changeList('')
fid1 = fopen(fileName,'r');
fid2 = fopen([fileName,'_new'],'w')
lineInFile = 1;
k = 1;
while lineInFile~=-1
  if ~isempty(strfind(lineInFile,'/seq/'))
    k = k+1;
  end
  fprintf(fid2,'%s\n',lineInFile);
  if k==5
    fprintf(fid2,'sleep 8h\n');
  end
  lineInFile = fgetl(fid1);
end

fclose(fid1);
fclose(fid2);

