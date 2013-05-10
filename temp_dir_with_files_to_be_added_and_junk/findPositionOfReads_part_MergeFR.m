function findPositionOfReads_part_MergeFR(dirName,sampleName,readLength,saveFileName,indexFirst,indexLast)
% rmpath(genpath('/home/csfaculty/shental/Weizmann/ph2users/fenoam/'));rmpath('/home/csfaculty/shental/matlab');cd ~/CS/mFilesBAC;addpath('~/CS/mFilesBAC');mcc -v -m findPositionOfReads_part_MergeFR.m -d /homes/csfaculty/shental/CS/BAC/cvx -a get_mean_posalign_Amnon_MergeForwardReverse.m

% findpositionofReads_part_MergeFR(auxData.dirNameForData,sampleName,readLength,saveFileName,1,20)

disp('run get_mean_posalign_MergeFR')

if ~isdeployed
  userDir = getuserdir;
else
  userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
  readLength = str2num(readLength);
  indexFirst = str2num(indexFirst);
  indexLast = str2num(indexLast);

end


load([userDir,'/CS/BAC/12Samples/Solexa/data/forAlign/primer750_primer_tail_charVector_afterIncludingRightTail']) 


data = load([userDir,'/CS/BAC/12Samples/Solexa/data/',dirName,'/data_',dirName,'_',sampleName]);
%keyboard

for i=indexFirst:indexLast
  
  % initialize
  out = get_mean_posalign_Amnon_MergeForwardReverse(data.uni_read_seq(i,:), database1,database_size);
  
  pos(i) = mean(out.I);
  sd(i) = std(out.I);
  rev_pos(i) = mean(bacseq_len(out.J)+1+1-out.I-readLength);
  rev_sd(i) = std(bacseq_len(out.J)+1+1-out.I-readLength);
  
end

save([saveFileName],'pos','rev_pos', 'sd', 'rev_sd','indexFirst','indexLast')
