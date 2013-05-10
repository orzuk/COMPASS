function runOneNodeForCompilation2(allocateToOneFirst,allocateToOneLast,basicAllocationFileName,saveFileName,userDir)

disp('runOneNodeForCompilation2.m The only difference is that it runs: prepareGroupOf1000DistributedSequenceFiles2 instead of prepareGroupOf1000DistributedSequenceFiles')

% the -a is important 
% rmpath(genpath('/home/csfaculty/shental/Weizmann/ph2users/fenoam/'));addpath('~/CS/mFilesBAC/');cd /homes/csfaculty/shental/CS/BAC/cvx;
% mcc -m runOneNodeForCompilation2 -a ./builtins -a ./commands -a ./functions -a ./keywords -a ./lib -a ./sdpt3 -a ./sedumi -a ./sets -a ./structures

if isdeployed
  allocateToOneFirst = str2num(allocateToOneFirst);
  allocateToOneLast = str2num(allocateToOneLast);
end
%addpath([userDir,'/CS/']);
%addpath(genpath([userDir,'/CS/BAC/cvx/']));
load(basicAllocationFileName)

for i=allocateToOneFirst:allocateToOneLast
  i
  [normalizedBac values] = prepareGroupOf1000DistributedSequenceFiles2(readLength,tmpInd{i},basicSeqNameDir,basicSeqKey);
  [fracRelevantReads,sumRelevantReads(i)] = currReads(uniqueReads,uniqueReads_length,values);clear values
  if sumRelevantReads(i)>0
    % find for 1000, output the sum of relevant reads 
    [x{i}]=runOneGroupOf1000ForCompilation(normalizedBac,fracRelevantReads,0);
  else
    x{i} = zeros(size(normalizedBac,2),1);
  end
end

save(saveFileName,'x','sumRelevantReads');
quit



