function prepareWindow_12Samples_createReads_MergeForwardAndReverse(auxData)
%keyboard
alignReads = load([auxData.uniqueReadsDirName,'/',auxData.uniqueReadsFilesName]); 
POS_data = load([auxData.POS_dirName,'/',auxData.POS_fileName,'_finalRes']);

readLengthPost = auxData.readLengthPost;
readLengthPre = auxData.readLengthPre;

uniqueReadsChar = alignReads.uni_read_seq;

%clear alignReads;

uni_reads = cell(size(uniqueReadsChar,1),1);
for i=1:size(uniqueReadsChar,1)
  uni_reads{i} = uniqueReadsChar(i,:);
end


r_post = repmat(char(zeros(size(uni_reads,1),readLengthPre)),readLengthPre-readLengthPost+1,1);
r_post = r_post(:,1:readLengthPost);

length_r_post_forward = zeros(size(r_post,1),1);
length_r_post_reverse = zeros(size(r_post,1),1);

r_POS = zeros(size(r_post,1),1);
r_REV_POS = zeros(size(r_post,1),1);
r_SD = zeros(size(r_post,1),1);
r_REV_SD = zeros(size(r_post,1),1);

r_shift = zeros(size(r_post,1),1);
k = 1;
for i=1:size(uni_reads,1)
  if mod(i,10^5)==1
    i
  end
  for j=1:readLengthPre-readLengthPost+1
    r_post(k,:) = uni_reads{i}(j:j+readLengthPost-1);
    
    length_r_post_forward(k) = alignReads.freq_ForwardAndReverse.forward(i);
    length_r_post_reverse(k) = alignReads.freq_ForwardAndReverse.reverse(i);

    r_shift(k) = j;
    
    
    r_POS(k) = POS_data.POS(i);
    r_REV_POS(k) = POS_data.REV_POS(i);
    r_SD(k) = POS_data.SD(i);
    r_REV_SD(k) = POS_data.REV_SD(i);
    
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

uniqueReads_length_forward = zeros(size(vals1,1),1);
uniqueReads_length_reverse = zeros(size(vals1,1),1);

shift_and_POS_data = cell(size(vals1,1),1);

for i=1:size(uniqueReads,1)
  uniqueReads_length_forward(i) = sum(length_r_post_forward(inds1{i}));
  uniqueReads_length_reverse(i) = sum(length_r_post_reverse(inds1{i}));
  
  shift_and_POS_data{i}.shift = r_shift(inds1{i});
  
  shift_and_POS_data{i}.POS = r_POS(inds1{i});
  shift_and_POS_data{i}.REV_POS = r_REV_POS(inds1{i});
  shift_and_POS_data{i}.SD = r_SD(inds1{i});
  shift_and_POS_data{i}.REV_SD = r_REV_SD(inds1{i});
    
  shift_and_POS_data{i}.length_forward = length_r_post_forward(inds1{i});
  shift_and_POS_data{i}.length_reverse = length_r_post_reverse(inds1{i});
  

  
end


save([auxData.saveDir,'/',auxData.fileName],'uniqueReads','shift_and_POS_data','uniqueReads_length_forward','uniqueReads_length_reverse')
