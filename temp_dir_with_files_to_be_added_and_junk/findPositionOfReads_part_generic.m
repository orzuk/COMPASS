function findPositionOfReads_part_generic(saveFileName,indexFirst,indexLast,charVectorFile,inputDirName,readsFileName,searchForwardAndReverseFlag)
% rmpath(genpath('/home/csfaculty/shental/Weizmann/ph2users/fenoam/'));rmpath('/home/csfaculty/shental/matlab');cd ~/CS/mFilesBAC;addpath('~/CS/mFilesBAC');mcc -v -m findPositionOfReads_part_generic.m -d /homes/csfaculty/shental/CS/BAC/cvx -a get_mean_posalign_generic 

if ~isdeployed
  userDir = getuserdir
else
  userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
  indexFirst = str2num(indexFirst);
  indexLast = str2num(indexLast);
  searchForwardAndReverseFlag = str2num(searchForwardAndReverseFlag);
end


load([charVectorFile]) 

load([inputDirName,'/',readsFileName]) 

for i=indexFirst:indexLast
  [pos(i) sd(i) rev_pos(i) rev_sd(i)] = get_mean_posalign_generic(uni_reads(i,:),database1,database_size,bacseq_len,searchForwardAndReverseFlag);
end

if searchForwardAndReverseFlag
  save([saveFileName],'pos','rev_pos', 'sd', 'rev_sd','indexFirst','indexLast')
else
  save([saveFileName],'pos','sd','indexFirst','indexLast')
end



 


