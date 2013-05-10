function prepareUniqueReads_Solexa(cellseq100_sample,readLength,sampleName,index)

for i=index
  clear red uniqueReads uniqueReads_length
  red = cellseq100_sample{i};
  
  found = [];
  for j=1:size(red,1)
    if ~isempty(find(red(j,:)~='A' & red(j,:)~='C'& red(j,:)~='G'& red(j,:)~='T'))
      found = [found,j];
      
    end
    if mod(j,10^5)==1
      disp(j)
    end
  end
  red(found,:) = [];
  
  red_packed = pack_seqs(red,64);
  clear red
  
  n = size(red_packed,1);
  groupSize = 5*10^5;
  part = 1:groupSize:n;
  part(end) = n+1;
  clear tmp_* junk*
  tmp_uniqueReads = [];
  tmp_uniqueReads_length = [];
  
  for j=1:length(part)-1
    tmp_red = red_packed(part(j):part(j+1)-1,:);
    [junk_uniqueReads,junk_uniqueReads_inds] = extract_sub_kmers(tmp_red, readLength*ones(size(tmp_red,1),1),readLength,1,0);
    [junk_junk_vals junk_junk_inds junk_uniqueReads_length] = get_duplicates(junk_uniqueReads_inds(:,1));
    tmp_uniqueReads = [tmp_uniqueReads;junk_uniqueReads];
    tmp_uniqueReads_length = [tmp_uniqueReads_length,junk_uniqueReads_length];
  end
  
  % unique of unique
  [uniqueReads,junk_uniqueReads_inds] = extract_sub_kmers(tmp_uniqueReads, readLength*ones(size(tmp_uniqueReads,1),1),readLength,1,0);
  [junk_junk_vals junk_junk_inds junk_junk_length] = get_duplicates(junk_uniqueReads_inds(:,1));
  
  clear red red_packed 
  
  uniqueReads_length = zeros(length(junk_junk_inds),1);
  for j=1:length(junk_junk_inds)
    uniqueReads_length(j) = sum(tmp_uniqueReads_length(junk_uniqueReads_inds(junk_junk_inds{j},2)));
  end
  clear tmp* junk*
  
  save(['~/CS/BAC/12Samples/Solexa/data/readsSolexa_WITH_REVERSE_',num2str(readLength),'_sample_',sampleName{i}],'uniqueReads','uniqueReads_length')
end




