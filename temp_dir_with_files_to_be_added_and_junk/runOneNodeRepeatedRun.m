function runOneNodeRepeatedRun(fileToLoad,l,mFileName)

fid = fopen([mFileName],'w');
while fid<0
  pause(5)
  fid = fopen([mFileName],'w');
end

fprintf(fid,'#!/bin/bash \n');
  
fprintf(fid,'TMPDIR=/tmp/${RANDOM}_${RANDOM}\n');
fprintf(fid,'mkdir $TMPDIR\n');
fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
fprintf(fid,'%s/CS/BAC/cvx/run_solveForGroupOr_Br.sh /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b %s %s %s %s %s %s\n',...
        userDir, ...
        num2str(allocateToOne(1)),' ',num2str(allocateToOne(end)) ,'a',basicAllocationFileName,saveFileName);
% '1' is just replacemnet for auxData needed in solveForGroup3                                    
fprintf(fid,'\n');
fprintf(fid,'rm -rf $TMPDIR \n');
fclose(fid);


for l=1:auxData.repeatRandomGroups
  [cX{l},cSumRelevantReads] = iterateParallelDistributedSequenceFilesOr(uniqueReads,uniqueReads_length,curr_kp,auxData.basicSeqNameDir,auxData.basicSeqKey,parallelType,[basicSaveName,'_k_',num2str(k)],userDir,[runFileName,'_k_',num2str(k)],auxData);
  n_kp(l,find(cX{l}>auxData.thresholdForCollectingBAC)) = 1;
  
  mx = max(mx,cX{l});
  mx_currSumRelevantReads = max(mx_currSumRelevantReads,cSumRelevantReads);
end

%%%%%%%%%%%%%%%%%%
load(fileToLoad)
l
[cX{l},cSumRelevantReads] = iterateParallelDistributedSequenceFilesOr(uniqueReads,uniqueReads_length,curr_kp,auxData.basicSeqNameDir,auxData.basicSeqKey,parallelType,[basicSaveName,'_k_',num2str(k),'_repeat_',num2str(l)],userDir,[runFileName,'_k_',num2str(k),'_repeat_',num2str(l)],auxData);


