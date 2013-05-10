function createUniqueReads_divideReads(fileNameInput,saveName,readLengthPre,readLengthPost)


load(fileNameInput)

r = int2nt(unpack_seqs(reads,readLengthPre,64));

readChar = cell(size(r,1),1);
for i=1:size(r,1)
  readChar{i} = r(i,:);
end

m4sort=sort(readChar);
[uni_reads,ia]=unique(m4sort,'first');
[~,ib]=unique(m4sort,'last');
basic_uniqueReads_length = ib-ia+1;

%keyboard

r_post = repmat(char(zeros(size(uni_reads,1),size(r,2))),readLengthPre-readLengthPost+1,1);
r_post = r_post(:,1:readLengthPost);
length_r_post = zeros(size(r_post,1),1);
k = 1;
for i=1:size(uni_reads,1)
  if mod(i,10^5)==1
    i
  end
  for j=1:readLengthPre-readLengthPost+1
    r_post(k,:) = uni_reads{i}(j:j+readLengthPost-1);
    length_r_post(k) = basic_uniqueReads_length(i);
    k = k+1;
  end
end
r_post(k+1:end,:) = [];


readChar = cell(size(r_post,1),1);
for i=1:size(r_post,1)
  readChar{i} = r_post(i,:);
end


[vals1 inds1 num_dups] = get_duplicates2(readChar);

uniqueReads = pack_seqs(vals1,64);
uniqueReads = cell2mat(uniqueReads);

uniqueReads_length = zeros(size(vals1,1),1);
for i=1:size(uniqueReads_length,1)
  uniqueReads_length(i) = sum(length_r_post(inds1{i}));
end

save(saveName,'uniqueReads_length','uniqueReads')



