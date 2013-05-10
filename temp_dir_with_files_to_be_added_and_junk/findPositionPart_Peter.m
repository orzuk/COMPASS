function findPositionPart_Peter(saveFileName,indexFirst,indexLast,dataPeterName)
% rmpath(genpath('/home/csfaculty/shental/Weizmann/ph2users/fenoam/'));rmpath('/home/csfaculty/shental/matlab');cd ~/CS/mFilesBAC;addpath('~/CS/mFilesBAC');mcc -v -m findPositionPart_Peter.m -d /homes/csfaculty/shental/CS/BAC/cvx -a get_mean_posalign_rna

if ~exist('dataPeterName')
  dataPeterName = 'dataPeter';
  disp('running with regular Peter in findPositionPart_Peter.m')
end


if ~isdeployed
  userDir = getuserdir
else
  userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
  indexFirst = str2num(indexFirst);
  indexLast = str2num(indexLast);

end


load([userDir,'/CS/BAC/12Samples/Solexa/data/forAlign_full16S/dataForAlign_full16S']) 

load([userDir,'/CS/BAC/Peter/',dataPeterName]) 
disp(['using ',dataPeterName])

for i=indexFirst:indexLast
  pos(i)=get_mean_posalign_rna(uni_reads{i},database1,database_size,bacseq_len);
end

save([saveFileName],'pos','indexFirst','indexLast')
