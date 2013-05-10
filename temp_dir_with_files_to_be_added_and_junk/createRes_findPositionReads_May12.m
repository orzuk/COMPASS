% with new get_mean_posalign_v2 - should compile into the generic one!!!!
% uses the only aligned reads


%%%%%%%%%
clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';

sampleNames = {'S1','S2','S3','S4','M1','M2','M3','M4','O7','O10','S7','S10'};

%%%%% read length 50
readLength = 100;
for i=[2:12]  
  sampleName = sampleNames{i};
  
  auxData = struct;
  auxData.queueName = 'hour';
  auxData.readLength = readLength;
  load([userDir,'/CS/BAC/12Samples/Solexa/data/opt3_th50_mar25_2012/illumina_reads',num2str(readLength),'_',sampleName,'_unifreq_corr_opt3_th50'],'freq_read')
  
  n = size(freq_read,1);
  clear freq_read
  
  auxData.numProcessors = round(n/750);
  auxData.numProcessors
  
  unix(['mkdir /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/findAlignment_May12/tmpRuns/alignTmp_',num2str(readLength),'_',sampleName])
  unix(['mkdir /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/findAlignment_May12/dataForSim/alignTmp_',num2str(readLength),'_',sampleName])
  
  findPositionOfReads_May12(sampleName,auxData);

end

