if 1==2 % creating reads

% sliding50 create reads
clear

dirName = 'sliding50';
createUniqueReads_12SamplesSolexa_sliding50(dirName,'O7',0)
createUniqueReads_12SamplesSolexa_sliding50(dirName,'O10',0)
createUniqueReads_12SamplesSolexa_sliding50(dirName,'S7',0)
createUniqueReads_12SamplesSolexa_sliding50(dirName,'S10',0)
createUniqueReads_12SamplesSolexa_sliding50(dirName,'M1',0)
createUniqueReads_12SamplesSolexa_sliding50(dirName,'M2',0)
createUniqueReads_12SamplesSolexa_sliding50(dirName,'M3',0)
createUniqueReads_12SamplesSolexa_sliding50(dirName,'M4',0)
createUniqueReads_12SamplesSolexa_sliding50(dirName,'S1',0)
createUniqueReads_12SamplesSolexa_sliding50(dirName,'S2',0)
createUniqueReads_12SamplesSolexa_sliding50(dirName,'S3',0)
createUniqueReads_12SamplesSolexa_sliding50(dirName,'S4',0)

%%%%%%%%%%%%%%%%%%55
  
end

if 1==2 % prepare the run files
  testData12SamplesSolexa_WITH_REVERSE_sliding50({'O7','O10','S7','S10','M1','M2','M3','M4','S1','S2','S3','S4'},'sliding50');
end

if 1==2
  createFilesFor_testData12Samples
end

if 1==2
  % mkdir resSolexa_WITH_REVERSE_sliding50
  % rm -f list_WITH_REVERSE_sliding50;find sol_WITH_REVERSE_sliding50* -name "*.mat" >> list_WITH_REVERSE_sliding50
  % basj
  % for i in `cat list_WITH_REVERSE_sliding50`; do cp $i resSolexa_WITH_REVERSE_sliding50;done
end

if 1==2
  sol_WITH_REVERSE_sliding50_S3
  sol_WITH_REVERSE_sliding50_O7
  sol_WITH_REVERSE_sliding50_S1
  sol_WITH_REVERSE_sliding50_S4
  % run again those that did not run
  testData12SamplesSolexa_WITH_REVERSE_sliding50_highCPU({'O7','S1','S3','S4'},'sliding50');
  
  createFilesFor_testData12Samples_highCPU
end

