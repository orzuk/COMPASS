function findPositionOfReads_part(sampleName,readLength,saveFileName,indexFirst,indexLast)
% rmpath(genpath('/home/csfaculty/shental/Weizmann/ph2users/fenoam/'));rmpath('/home/csfaculty/shental/matlab');cd ~/CS/mFilesBAC;addpath('~/CS/mFilesBAC');mcc -v -m findPositionOfReads_part.m -d /homes/csfaculty/shental/CS/BAC/cvx -a get_mean_posalign_v2


if ~isdeployed
  userDir = getuserdir
else
  userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
  readLength = str2num(readLength);
  indexFirst = str2num(indexFirst);
  indexLast = str2num(indexLast);

end


load([userDir,'/CS/BAC/12Samples/Solexa/data/forAlign/primer750_primer_tail_charVector_afterIncludingRightTail']) 

load([userDir,'/CS/BAC/12Samples/Solexa/data/primer_tail/illumina_reads',num2str(readLength),'_',sampleName,'_uni_primer_tail_noR']) 

w = ['uni_reads50 = uni_reads',num2str(readLength),'_',sampleName,';'];
eval(w);



for i=indexFirst:indexLast
  pos(i)=get_mean_posalign_v2(uni_reads50(i,:),database1,database_size,bacseq_len);
end

save([saveFileName],'pos','indexFirst','indexLast')
