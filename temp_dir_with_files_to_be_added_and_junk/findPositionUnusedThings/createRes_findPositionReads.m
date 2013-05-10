
clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';

sampleNames = {'O7','O10','S7','S10','M1','M2','M3','M4','S1','S2','S3','S4'};
% 11,5
for i=[1:4 6:10 12]
  sampleName = sampleNames{i};
  
  auxData = struct;auxData.queueName = 'hour';
  load([userDir,'/CS/BAC/12Samples/Solexa/data/primer_tail/illumina_reads50_',sampleName,'_uni'],'freq*')
  w = ['freq = freq_uni_reads50_',sampleName,';'];
  eval(w);
  n = size(freq,1);
  
  auxData.numProcessors = round(n/2000);
  auxData.numProcessors
  
  unix(['mkdir /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/findAlignment/tmpRuns/alignTmp_',sampleName])
  unix(['mkdir /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/findAlignment/dataForSim/alignTmp_',sampleName])
  
  findPositionOfReads(sampleName,auxData);

end

% do S7 again (3)

% 1 2 5 11

% left: 4 6 7 8 9 10 12
for i=[4 6 7]
  
for i=[8 9 10]
  
for i=[10 12]
  



%%%%%%%%%%%%%55
% run in batch

sampleNames = {'O7','O10','S7','S10','M1','M2','M3','M4','S1','S2','S3','S4'};
% 11,5

