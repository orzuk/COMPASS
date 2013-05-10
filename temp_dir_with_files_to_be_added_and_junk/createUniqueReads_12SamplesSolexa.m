function createUniqueReads_12SamplesSolexa(dirName,samplesName,readLength,addOneFlag)

if ~exist('addOneFlag')
  addOneFlag = 0;
end
load(['~/CS/BAC/12Samples/Solexa/data/',dirName,'/illumina_reads',num2str(readLength),'_',samplesName,'_uni']);


w = ['uniqueReads_length = freq_uni_reads',num2str(readLength),'_',samplesName,';'];
eval(w);

if addOneFlag
  disp('added 1. createUniqueReads_12SamplesSolexa.m')
  uniqueReads_length = uniqueReads_length + 1;
else
  disp('does not add 1. createUniqueReads_12SamplesSolexa.m')
end

uniqueReads_length = uniqueReads_length';

w = ['tmpReads = cell2mat(uni_reads',num2str(readLength),'_',samplesName,');'];
eval(w)

uniqueReads = pack_seqs(tmpReads,64);

save(['~/CS/BAC/12Samples/Solexa/data/',dirName,'/readsSolexa_WITH_REVERSE_',num2str(readLength),'_sample_',samplesName],'uniqueReads','uniqueReads_length');


%keyboard