function [fam,val]=readClassifedReads(fileName,num_unique_reads,typIndex)

if ~exist('typIndex')
  typIndex = 'fromFile'
else
  typIndex = 'roling';
end

%keyboard
fam = cell(num_unique_reads,1);
val = zeros(num_unique_reads,1);


fid = fopen(fileName,'r');
lineInFile = ' ';
while isempty(strfind(lineInFile,'Details'))
    lineInFile = fgetl(fid);
end

lineInFile = fgetl(fid);

if strcmp(typIndex,'fromFile')
  while lineInFile~=-1 
    a = findstr(lineInFile,';');
    num = str2num(lineInFile(1:a(1)-1));
    genus = lineInFile(a(6)+1:a(7)-1);
    fam{num} = genus;
    b = lineInFile(a(7)+1:a(8)-1);
    c = find(b=='%');
    val(num) = str2num(b(1:c-1));
    lineInFile = fgetl(fid);
    %num
  end
  
else
  num = 1;
  while lineInFile~=-1 
    a = findstr(lineInFile,';');
    genus = lineInFile(a(6)+1:a(7)-1);
    fam{num} = genus;
    b = lineInFile(a(7)+1:a(8)-1);
    c = find(b=='%');
    val(num) = str2num(b(1:c-1));
    lineInFile = fgetl(fid);
    num = num+1;
    %num
  end
end
fclose(fid); 



