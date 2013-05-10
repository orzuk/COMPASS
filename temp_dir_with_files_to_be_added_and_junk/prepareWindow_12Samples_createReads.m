function prepareWindow_12Samples_createReads(auxData)
%keyboard
alignReads = load([auxData.inputDirName,'/',auxData.uniqueReadsFilesName],'uni_read_seq','freq_read','POS'); 
readLengthPost = auxData.readLengthPost;
readLengthPre = auxData.readLengthPre;

uniqueReadsChar = alignReads.uni_read_seq;
uniqueReads_length = alignReads.freq_read;
POS = alignReads.POS;

clear alignReads;

uni_reads = cell(size(uniqueReadsChar,1),1);
for i=1:size(uniqueReadsChar,1)
  uni_reads{i} = uniqueReadsChar(i,:);
end


r_post = repmat(char(zeros(size(uni_reads,1),readLengthPre)),readLengthPre-readLengthPost+1,1);
r_post = r_post(:,1:readLengthPost);
length_r_post = zeros(size(r_post,1),1);
r_POS = zeros(size(r_post,1),1);
r_shift = zeros(size(r_post,1),1);
k = 1;
for i=1:size(uni_reads,1)
  if mod(i,10^5)==1
    i
  end
  for j=1:readLengthPre-readLengthPost+1
    r_post(k,:) = uni_reads{i}(j:j+readLengthPost-1);
    length_r_post(k) = uniqueReads_length(i);
    r_shift(k) = j;
    r_POS(k) = POS(i);
    k = k+1;
  end
end


uni_reads = cell(size(r_post,1),1);
for i=1:size(r_post,1)
  uni_reads{i} = r_post(i,:);
end


[vals1 inds1 num_dups] = get_duplicates2(uni_reads);

uniqueReads = pack_seqs(vals1,64);
uniqueReads = cell2mat(uniqueReads);

uniqueReads_length = zeros(size(vals1,1),1);
shift_and_POS_data = cell(size(vals1,1),1);

for i=1:size(uniqueReads_length,1)
  uniqueReads_length(i) = sum(length_r_post(inds1{i}));
  shift_and_POS_data{i}(1,:) = r_shift(inds1{i})';
  shift_and_POS_data{i}(2,:) = r_POS(inds1{i})';
  shift_and_POS_data{i}(3,:) = length_r_post(inds1{i});
end


save([auxData.saveDir,'/',auxData.fileName],'uniqueReads','uniqueReads_length','shift_and_POS_data')
