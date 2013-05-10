function createUniqueReads_12SamplesSolexa_sliding50(dirName,samplesName,addOneFlag)

if ~exist('addOneFlag')
  addOneFlag = 0;
end
load(['~/CS/BAC/12Samples/Solexa/data/',dirName,'/illumina_reads50slide100','_',samplesName,'_uni']);

%keyboard
w = ['uniqueReads_length = reads_uni_freq;'];
eval(w);

if addOneFlag
  disp('added 1. createUniqueReads_12SamplesSolexa.m')
  uniqueReads_length = uniqueReads_length + 1;
else
  disp('does not add 1. createUniqueReads_12SamplesSolexa_sliding50.m')
end

uniqueReads_length = uniqueReads_length';
clear reads_uni_freq

uniqueReads = pack_seqs(reads_uni,64);

save(['~/CS/BAC/12Samples/Solexa/data/',dirName,'/readsSolexa_WITH_REVERSE_sliding50_sample_',samplesName],'uniqueReads','uniqueReads_length');


%keyboard