function createUniqueReads_12SamplesSolexa_eq_unicov_MergeForwardReverse(dirName,readLength,samplesName,typeName,thresholdForReverse)

load(['~/CS/BAC/12Samples/Solexa/data/',dirName,'/illumina_reads',num2str(readLength),'_',samplesName,typeName],'uni_read_seq','freq_read','freq_corr2pos','POS');

disp(['thresholdForReverse is taken as: ',num2str(thresholdForReverse)])


% see if there are non ACGT
if 1==2
  n = zeros(size(uni_read_seq,1),1);
  for i=1:size(uni_read_seq,1)
    if ~isempty(find(uni_read_seq(i,:)~='A' & uni_read_seq(i,:)~='C' & uni_read_seq(i,:)~='G' & uni_read_seq(i,:)~='T'))
      n(i) = 1;
    end
  end
  length(find(n))
  return
end

%%%%%%%%%%%%

if length(unique(uni_read_seq))~=4
  disp('non ACGT. createUniqueReads_12SamplesSolexa_eq_unicov_MergeForwardReverse.m')
  keyboard
end

%keyboard

readChar = mat2cell(uni_read_seq,ones(size(uni_read_seq,1),1),size(uni_read_seq,2));
[vals1 inds1 num_dups] = get_duplicates2(readChar);

if max(num_dups)~=2
  disp('problem. createUniqueReads_12SamplesSolexa_eq_unicov_MergeForwardReverse.m')
  keyboard
end


new_uni_reads = cell2mat(vals1);

% reads again
sum_freq = zeros(1,length(num_dups));
freq_ForwardAndReverse = struct;
freq_ForwardAndReverse.forward = zeros(1,length(num_dups));
freq_ForwardAndReverse.reverse = zeros(1,length(num_dups));
for i=1:length(num_dups)
  sum_freq(i) = sum(freq_corr2pos(inds1{i}));
  for j=1:length(inds1{i})
    if POS(inds1{i}(j))>thresholdForReverse
      freq_ForwardAndReverse.forward(i) = freq_read(inds1{i}(j));
    else
      freq_ForwardAndReverse.reverse(i) = freq_read(inds1{i}(j));
    end
  end
end

uniqueReads = pack_seqs(new_uni_reads,64);

% with equalization
uniqueReads_length = sum_freq;
uni_read_seq = new_uni_reads;
save(['~/CS/BAC/12Samples/Solexa/data/',dirName,'/data_',dirName,'_',samplesName],'uniqueReads','uniqueReads_length','uni_read_seq','freq_ForwardAndReverse');
%save(['~/CS/BAC/12Samples/Solexa/data/',dirName,'/data_',dirName,'_',samplesName],'uniqueReads','uni_read_seq','freq_ForwardAndReverse');

%keyboard
