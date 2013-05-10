function createUniqueReads_12SamplesSolexa_with_or_without_equalization(dirName,readLength,samplesName)

load(['~/CS/BAC/12Samples/Solexa/data/',dirName,'/illumina_reads',num2str(readLength),'_',samplesName,'_uni_corrfreq'],'uni_read_seq','freq_read','freq_corr2pos');

%keyboard

uniqueReads = pack_seqs(uni_read_seq,64);

% no equalization
uniqueReads_length = freq_read';
save(['~/CS/BAC/12Samples/Solexa/data/',dirName,'/readsSolexa_without_equalization_',num2str(readLength),'_',samplesName],'uniqueReads','uniqueReads_length');

% with equalization
uniqueReads_length = freq_corr2pos';
save(['~/CS/BAC/12Samples/Solexa/data/',dirName,'/readsSolexa_with_equalization_',num2str(readLength),'_',samplesName],'uniqueReads','uniqueReads_length');

%keyboard