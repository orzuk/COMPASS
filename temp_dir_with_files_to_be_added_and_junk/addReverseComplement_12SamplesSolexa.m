function addReverseComplement_12SamplesSolexa(sampleName,readLength)
%addReverseComplement_12SamplesSolexa({'O7','O10','S7','S10','M1','M2','M3','M4','S1','S2','S4'},100);addReverseComplement_12SamplesSolexa({'O7','O10','S7','S10','M1','M2','M3','M4','S1','S2','S4'},50)

for i=1:length(sampleName)
  %keyboard
  clear uniqueReads uniqueReads_length
  load(['~/CS/BAC/12Samples/Solexa/data/readsSolexa_',num2str(readLength),'_sample_',sampleName{i}]);
  
  seq = unpack_seqs(uniqueReads,readLength,64);
  seqRev = pack_seqs(seqrcomplement(seq),64);;
  %keyboard
  uniqueReads = [uniqueReads;seqRev];
  uniqueReads_length = [uniqueReads_length;uniqueReads_length];
    
  save(['~/CS/BAC/12Samples/Solexa/data/readsSolexa_WITH_REVERSE_',num2str(readLength),'_sample_',sampleName{i}],'uniqueReads','uniqueReads_length')
    
  sampleName{i},size(uniqueReads_length)
end  

