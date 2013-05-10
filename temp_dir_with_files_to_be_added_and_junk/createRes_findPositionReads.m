%%%%%%%
% prepare the data
clear
load ~/CS/BAC/12Samples/Solexa/data/forAlign/bac16s_primers750_primer_tail_chararray

%Sequence_uni_amp_chararray = Sequence_uni_amp_chararray';
database1 = Sequence_uni_amp_chararray(:);
database_size = size(Sequence_uni_amp_chararray);




save ~/CS/BAC/12Samples/Solexa/data/forAlign/primer750_primer_tail_charVector_afterIncludingRightTail database1 database_size bacseq_len


%%%%%%%%%%%%


clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';

sampleNames = {'O7','O10','S7','S10','M1','M2','M3','M4','S1','S2','S3','S4'};

%%%%% read length 50
readLength = 50;
for i=[12:-1:1]  %11
  sampleName = sampleNames{i};
  
  auxData = struct;
  auxData.queueName = 'hour';
  auxData.readLength = readLength;
  load([userDir,'/CS/BAC/12Samples/Solexa/data/primer_tail/illumina_reads',num2str(readLength),'_',sampleName,'_uni_primer_tail_noR'],'freq*')
  w = ['freq = freq_uni_reads',num2str(readLength),'_',sampleName,';'];
  eval(w);
  n = size(freq,1);
  
  auxData.numProcessors = round(n/1500);
  auxData.numProcessors
  
  unix(['mkdir /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/findAlignment/tmpRuns/alignTmp_',num2str(readLength),'_',sampleName])
  unix(['mkdir /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/findAlignment/dataForSim/alignTmp_',num2str(readLength),'_',sampleName])
  
  findPositionOfReads(sampleName,auxData);

end


%%%%%%%%%
clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';

sampleNames = {'O7','O10','S7','S10','M1','M2','M3','M4','S1','S2','S3','S4'};

%%%%% read length 50
readLength = 100;
for i=[12] 
  sampleName = sampleNames{i};
  
  auxData = struct;
  auxData.queueName = 'hour';
  auxData.readLength = readLength;
  load([userDir,'/CS/BAC/12Samples/Solexa/data/primer_tail/illumina_reads',num2str(readLength),'_',sampleName,'_uni_primer_tail_noR'],'freq*')
  w = ['freq = freq_uni_reads',num2str(readLength),'_',sampleName,';'];
  eval(w);
  n = size(freq,1);
  
  auxData.numProcessors = round(n/1500);
  auxData.numProcessors
  
  unix(['mkdir /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/findAlignment/tmpRuns/alignTmp_',num2str(readLength),'_',sampleName])
  unix(['mkdir /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/findAlignment/dataForSim/alignTmp_',num2str(readLength),'_',sampleName])
  
  findPositionOfReads(sampleName,auxData);

end


%%%%%%%%%%%%%55
% run in batch

sampleNames = {'O7','O10','S7','S10','M1','M2','M3','M4','S1','S2','S3','S4'};
% 11,5

