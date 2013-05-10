% check the reverse issue
if 1==2
  clear
  sampleNames = {'S1','S2','S3','S4','M1','M2','M3','M4','O7','O10','S7','S10'};
  
  i=1
  dirNameForData = 'opt3_MergeForwardReverse';
  u = load(['~/CS/BAC/12Samples/Solexa/data/',dirNameForData,'/data_',dirNameForData,'_',sampleNames{i}]);
  r = load(['~/CS/BAC/12Samples/Solexa/data/primer_tail/illumina_reads100_S1_uni_primer_tail_noR'])
  

  seq = u.uni_read_seq(1,:);
  
  seqRev = seqrcomplement(seq);
  
  l_f = 0;
  l_r = 0;
  for j=1:size(r.freq_uni_reads100_S1)
    if strcmp(r.uni_reads100_S1(j,:),seq);
      l_f = l_f+r.freq_uni_reads100_S1(j);
    end
    if strcmp(r.uni_reads100_S1(j,:),seqRev);
      l_r = l_r+r.freq_uni_reads100_S1(j);
    end
  end
  
  % if marked reverse than appears are reverse in the original reads
end



%%%%%%%%%%%%
% align the MergeForwardAndReverse files

%%%%%%%%%
clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';

sampleNames = {'S1','S2','S3','S4','M1','M2','M3','M4','O7','O10','S7','S10'};
dirNameForData = 'opt3_MergeForwardReverse';
%%%%% read length 50
readLength = 100;
for i=11:12%[1:12]  
  sampleName = sampleNames{i};
  
  auxData = struct;
  auxData.queueName = 'hour';
  auxData.readLength = readLength;  
  auxData.dirNameForData = dirNameForData;
  u = load([userDir,'/CS/BAC/12Samples/Solexa/data/',dirNameForData,'/data_',dirNameForData,'_',sampleNames{i}],'uniqueReads_length');

  
  n = size(u.uniqueReads_length,2);
  clear u
  
  auxData.numProcessors = round(n/500);
  auxData.numProcessors 
  
  unix(['mkdir /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/findAlignment_MergeForwardAndReverse/tmpRuns/alignTmp_',num2str(readLength),'_',sampleName])
  unix(['mkdir /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/findAlignment_MergeForwardAndReverse/dataForSim/alignTmp_',num2str(readLength),'_',sampleName])
  
  
  findPositionOfReads_MergeForwardAndReverse(sampleName,auxData)
end

%%%%%%%%%%%%%%%%%%%%%%%%%555555
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end prepare the positions for reads

% send Amnon the reads and their position

%%%%%%%%%%%%%%%%%%%%%%%%55
% run the normalized case


