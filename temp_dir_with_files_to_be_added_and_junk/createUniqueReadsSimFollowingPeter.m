function createUniqueReadsSimFollowingPeter(fileNameInput,saveName,readLength)


load(fileNameInput)
%keyboard
r = int2nt(unpack_seqs(reads,readLength,64));
readChar = cell(size(r,1),1);
for i=1:size(r,1)
  readChar{i} = r(i,:);
end

m4sort=sort(readChar);
[uni_reads,ia]=unique(m4sort,'first');
[~,ib]=unique(m4sort,'last');
uniqueReads_length = ib-ia+1;

uniqueReads = pack_seqs(uni_reads,64);
uniqueReads = cell2mat(uniqueReads);

save(saveName,'uniqueReads_length','uniqueReads')
