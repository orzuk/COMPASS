function createRes_findPositionReads_batch(sampleName)
% rmpath(genpath('/home/csfaculty/shental/Weizmann/ph2users/fenoam/'));rmpath('/home/csfaculty/shental/matlab');cd ~/CS/mFilesBAC;addpath('~/CS/mFilesBAC');mcc -v -m createRes_findPositionReads_batch.m -d /homes/csfaculty/shental/CS/BAC/cvx -a findPositionOfReads

userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';

auxData = struct;auxData.queueName = 'hour';
load([userDir,'/CS/BAC/12Samples/Solexa/data/primer_tail/illumina_reads50_',sampleName,'_uni'],'freq*')
w = ['freq = freq_uni_reads50_',sampleName,';'];
eval(w);
n = size(freq,1);

auxData.numProcessors = round(n/2000);
auxData.numProcessors
  
findPositionOfReads(sampleName,auxData);




