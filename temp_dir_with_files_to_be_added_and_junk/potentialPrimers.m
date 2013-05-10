%%%%%%%%%%%%%%%%%
% find frequent 20 mers 
clear
userDir = getuserdir;
addpath([userDir,'/CS/mFilesBAC/']);
addpath(genpath([userDir,'/CS/BAC/cvx/']))

basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];
basicSeqKey= [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];

curr_kp = 1:410849;

n = length(curr_kp);
part = 1:1000:n;
part(end) = n+1;

tmpInd = cell(length(part)-1,1);
for i=1:length(part)-1
  tmpInd{i} = part(i):part(i+1)-1;
  tmpInd{i} = curr_kp((tmpInd{i}));
end

readLength = 20;

pot = cell(1,length(tmpInd));

matlabpool open local 6
parfor i=1:length(tmpInd)
  i
  [normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOrFourth(readLength,tmpInd{i},basicSeqNameDir,basicSeqKey); 
  
  a = sum(normalizedBac,2);
  b = find(a>0.5*size(normalizedBac,2) & a>prctile(a,95));

  pot{i} = [values(b,:) full(a(b))];
end
matlabpool close


vals = cell2mat(pot');

[junk_vals junk_inds uniqueReads_length] = get_duplicates(vals(:,1));


list(:,1) = junk_vals;
for i=1:size(list,1)
  list(i,2) = sum(vals(junk_inds{i},2));
end

numOccurances = list(:,2);
seqs = int2nt(unpack_seqs(list(:,1), readLength, 64));

save ~/CS/BAC/occurancesOf20Mers numOccurances seqs
