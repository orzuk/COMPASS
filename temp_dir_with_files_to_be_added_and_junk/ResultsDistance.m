% calculate the distance matrix for a reconstruction solution file and save
% it
% input:
% dirName,fileName - the input directory and filename
% output:
% a saved file (like the input file but with a 'distmat-' beginning, which
% contains the following matrices (of size orig * rec):
% idenFracDist - fraction of mismatches (mismatches/minimum length)
% idenFracDistMax - fraction of mismatches (mismatches/maximum length)
% scoreDist - the pairwise alignment scores

function [idenFracDist,idenFracDistMax,scoreDist]=ResultsDistance(dirName,fileName,auxData)

% rmpath(genpath('/home/csfaculty/shental/Weizmann/ph2users/fenoam/'));addpath('/home/csfaculty/shental/CS/mFilesBAC/');cd /home/csfaculty/shental/CS/BAC/cvx; mcc -m ResultsDistance -a MahDist



if ~isfield(auxData,'brFlag')
  auxData.brFlag = 0;
end
if auxData.brFlag
  userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
  disp('run in Br')
  parallelType = 'cluster';
else
  userDir = getuserdir;
  [junk systemName] = system('hostname');
  systemName = deblank(systemName);
  if ~isempty(findstr(userDir,'fenoam')) & ~strcmp(systemName,'complex1') %WIS
    parallelType = 'cluster';
    disp('runs in WIS')
  elseif ~isempty(findstr(userDir,'fenoam')) & strcmp(systemName,'complex1') %WIS
    parallelType = 'local';
    disp('run in complex1')
  else %OU
    parallelType = 'local';
    disp('run in OU')
  end
end

if ~isfield(auxData,'basicSeqNameDir') 
  auxData.basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];
end
if ~isfield(auxData,'auxData.basicSeqKey')
  auxData.basicSeqKey = [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
end
disp(['working with database: ',auxData.basicSeqNameDir,'; Is it ok?'])



threshold=0;

r = load([dirName fileName]);

orig=r.resCell{1};
rec=r.resCell{2};

origSeqs=orig(:,1);
recSeqs=rec(:,1);
recFreq=rec(:,2);

[idenFracDist,idenFracDistMax,scoreDist]=MahDist(origSeqs,recSeqs,recFreq,threshold,auxData);

save([dirName 'distmat_' fileName],'idenFracDist','idenFracDistMax','scoreDist','orig','rec');
