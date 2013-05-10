function runOneNodeForCompilation4(allocateToOneFirst,allocateToOneLast,basicAllocationFileName,saveFileName,userDir)

disp('runOneNodeForCompilation2.m The only difference is that it runs: prepareGroupOf1000DistributedSequenceFiles2 instead of prepareGroupOf1000DistributedSequenceFiles')

% the -a is important 
% rmpath(genpath('/home/csfaculty/shental/Weizmann/ph2users/fenoam/'));addpath('~/CS/mFilesBAC/');cd /homes/csfaculty/shental/CS/BAC/cvx;
% mcc -m runOneNodeForCompilation4 -a ./builtins -a ./commands -a ./functions -a ./keywords -a ./lib -a ./sdpt3 -a ./sedumi -a ./sets -a ./structures

if isdeployed
  allocateToOneFirst = str2num(allocateToOneFirst);
  allocateToOneLast = str2num(allocateToOneLast);
end
%addpath([userDir,'/CS/']);
%addpath(genpath([userDir,'/CS/BAC/cvx/']));

[x,sumRelevantReads,independentInd,dependentInd] = solveForGroup4(allocateToOneFirst,allocateToOneLast,auxData,basicAllocationFileName);

save(saveFileName,'x','sumRelevantReads','independentInd','dependentInd');
quit



