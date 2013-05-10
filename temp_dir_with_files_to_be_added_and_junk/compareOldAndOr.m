
% reads

% red are the reads

[uniqueReads,uniqueReads_inds] = extract_sub_kmers(red, auxData.readLength*ones(size(red,1),1),auxData.readLength, 1,0);
  
[vals inds num_dups] = get_duplicates(uniqueReads_inds(:,1));
  
uniqueReads0 = 'n'*zeros(size(uniqueReads,1),50,'uint8');
for l=1:size(uniqueReads,1)
  uniqueReads0(l,:) = int2nt(unpack_seqs(uniqueReads(l,:),50,64));
end
   
  
  
% check that they are the same
red1 = zeros(size(red,1),50);
for l=1:size(red,1)
  red1(l,:) = int2nt(unpack_seqs(red(l,:),50,64));
end
[uniqueReads1,uniqueReads_inds1] = myUniqueUINT8(red1);
 

u0 = sortrows(uniqueReads0);
u1 = sortrows(uniqueReads1);
  
find(u0-u1)


%%%%%%%%%%%%%%%
% in prepareGroupOf1000DistributedSequenceFilesOr

[normalizedBac values] = BuildMixingMatrixFromSequences(readLength,Sequence1,len_uni(tmpInd));

[normalizedBac1 values1] = prepareGroupOf1000DistributedSequenceFiles2(readLength,tmpInd,'/home/csfaculty/shental/CS/BAC/datNoNonACGT/unpacked/',basicSeqKey);

size(normalizedBac)
size(normalizedBac1)

for j=1:length(Sequence1)
  a = int2nt(unpack_seqs(Sequence1{j},len_uni(tmpInd(j)),64));
  if ~isempty(find(a~='A' & a~='C'  & a~='G'  & a~='T' ))
    a
  end
end

%%%%%%%%%%55
% is solveForGroupOr - compare currReads - finite
[fracRelevantReads,sumRelevantReads(i)] = currReads(uniqueReads,uniqueReads_length,values,auxData.inifiniteNumberOfReadsFlag,dataIn);

uniqueReads1 = zeros(size(uniqueReads,1),50);
for l=1:size(uniqueReads,1)
  uniqueReads1(l,:) = int2nt(unpack_seqs(uniqueReads(l,:),50,64));
end
values1 = zeros(size(values,1),50,'uint8');
for l=1:size(values,1)
  if mod(l,1000)==1
    l
  end
  values1(l,:) = int2nt(unpack_seqs(values(l,:),50,64));
end

[fracRelevantReads1,sumRelevantReads1] = currReads(uniqueReads1,uniqueReads_length,values1,auxData.inifiniteNumberOfReadsFlag,dataIn);

find(fracRelevantReads1-fracRelevantReads)


%%%%%%%%5
% is solveForGroupOr - compare currReads - infinite

