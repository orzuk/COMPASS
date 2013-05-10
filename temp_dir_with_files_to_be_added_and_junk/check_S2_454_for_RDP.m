clear
load /homes/csfaculty/shental/CS/BAC/12Samples/454/data/reads_samples_S2


reads = int2nt(unpack_seqs(uniqueReads,350,64));

Header = cell(1,size(uniqueReads,1));
for i=1:length(Header)
  Header{i} = [num2str(i),' unique: ',num2str(uniqueReads_length(i))];
end

fastawrite('~/CS/BAC/12Samples/454/data/check_S2/data_reads_350_S2.fa',Header,reads);

fid = fopen('~/CS/BAC/12Samples/454/data/check_S2/allrank_data_reads_350_S2.fa_classified.txt','r');
lineInFile = ' ';

while isempty(strfind(lineInFile,'Symbol +/-'))
    lineInFile = fgetl(fid);
end
lineInFile = fgetl(fid);
lineInFile = fgetl(fid);

count = 0;
while lineInFile~=-1 
   lineInFile = fgetl(fid);
   a = find(lineInFile==';');
   num = str2num(lineInFile(1:a(1)-1));
   if ~isempty(findstr(lineInFile,'Lacto'))
     lineInFile
     %pause
     count = count+uniqueReads_length(num);
   end
end
fclose(fid)

count/sum(uniqueReads_length)


%%%%%%%%%%%%%%%%%%%
% do the same for S1

clear
fclose('all')
load /homes/csfaculty/shental/CS/BAC/12Samples/454/data/reads_samples_S1


reads = int2nt(unpack_seqs(uniqueReads,350,64));

Header = cell(1,size(uniqueReads,1));
for i=1:length(Header)
  Header{i} = [num2str(i),' unique: ',num2str(uniqueReads_length(i))];
end

fastawrite('~/CS/BAC/12Samples/454/data/check_S1/data_reads_350_S1.fa',Header,reads);

clear lineInFile
fid = fopen('~/CS/BAC/12Samples/454/data/check_S1/allrank_data_reads_350_S1.fa_classified.txt','r');
lineInFile = ' ';

while isempty(strfind(lineInFile,'Symbol +/-'))
    lineInFile = fgetl(fid);
end
lineInFile = fgetl(fid);
lineInFile = fgetl(fid);

count = 0;
while lineInFile~=-1 
   a = find(lineInFile==';');
   num = str2num(lineInFile(1:a(1)-1));
   lineInFile
   if ~isempty(findstr(lineInFile,'Lacto'))
     lineInFile
     %pause
     count = count+uniqueReads_length(num);
   end
   lineInFile = fgetl(fid);
end
fclose(fid)


