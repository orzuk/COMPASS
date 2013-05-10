if 1==2 % creating reads

% create reads
clear

dirName = 'with_eq_unicov_noam_mar21';
typeName = '_unicov_corrfreq_noam';

readLength = 100;
createUniqueReads_12SamplesSolexa_with_eq_unicov_corrfreq(dirName,readLength,'O7',typeName)
createUniqueReads_12SamplesSolexa_with_eq_unicov_corrfreq(dirName,readLength,'O10',typeName)
createUniqueReads_12SamplesSolexa_with_eq_unicov_corrfreq(dirName,readLength,'S7',typeName)
createUniqueReads_12SamplesSolexa_with_eq_unicov_corrfreq(dirName,readLength,'S10',typeName)
createUniqueReads_12SamplesSolexa_with_eq_unicov_corrfreq(dirName,readLength,'M1',typeName)
createUniqueReads_12SamplesSolexa_with_eq_unicov_corrfreq(dirName,readLength,'M2',typeName)
createUniqueReads_12SamplesSolexa_with_eq_unicov_corrfreq(dirName,readLength,'M3',typeName)
createUniqueReads_12SamplesSolexa_with_eq_unicov_corrfreq(dirName,readLength,'M4',typeName)
createUniqueReads_12SamplesSolexa_with_eq_unicov_corrfreq(dirName,readLength,'S1',typeName)
createUniqueReads_12SamplesSolexa_with_eq_unicov_corrfreq(dirName,readLength,'S2',typeName)
createUniqueReads_12SamplesSolexa_with_eq_unicov_corrfreq(dirName,readLength,'S3',typeName)
createUniqueReads_12SamplesSolexa_with_eq_unicov_corrfreq(dirName,readLength,'S4',typeName)

% move with and without to their directories


%%%%%%%%%%%%%%%%%%55
  
end

if 1==2 % prepare the run files equalization
  
  % with equalization new format
  testData12Samples_with_eq_unicov_corrfreq({'S1','O7','O10','S7','S10','M1','M2','M3','M4','S2','S3','S4'},'with_eq_unicov_noam_mar21','with_eq_unicov_noam_mar21')
end

if 1==2
  % create files for run
  createFilesFor_testData12Samples_with_eq_unicov('with_eq_unicov_noam_mar21')
end


if 1==2
  mkdir resSolexa_with_eq_unicov_noam_mar21
  % rm -f list_with_eq_unicov_noam_mar21;find sol_with* -name "*.mat" >> list_with_eq_unicov_noam_mar21
  % basj
  % for i in `cat list_with_eq_unicov_noam_mar21`; do cp $i resSolexa_with_eq_unicov_noam_mar21;done
end


if 1==2
  % create the matrix 
  clear
  load ~/CS/BAC/resSolexa_with_eq_unicov_mar21/sol_with_eq_unicov_noam_mar21_S1_noCorrection_100
  load /home/csfaculty/shental/CS/BAC/12Samples/Solexa/data/with_eq_unicov_noam_mar21/readsSolexa_with_equalization_with_eq_unicov_noam_mar21_100_S1
  readLength = 100;
  
  userDir = getuserdir;
  datasetName = 'primers750_primer_tail';
  basicSeqNameDir = [userDir,'/CS/BAC/',datasetName,'/datNoNonACGT/packed64/'];
  basicSeqKey = [userDir,'/CS/BAC/',datasetName,'/datNoNonACGT/keyNoNonACGT_',datasetName];
  
   
  tmpInd = find(abs(found{3})>10^-2);
  [normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOrFourth(readLength,tmpInd,basicSeqNameDir,basicSeqKey);
  
  
  dataIn = struct;
  [fracRelevantReads,sumRelevantReads] = currReadsFourth(uniqueReads,uniqueReads_length,values,0,dataIn);
  
  
  
  save ~/CS/BAC/dataForS1 fracRelevantReads found tmpInd normalizedBac values uniqueReads uniqueReads_length
  

end