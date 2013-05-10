if 1==2 % creating reads

% sliding50 create reads
clear

dirName = 'with_or_without_equalization';
readLength = 50;
createUniqueReads_12SamplesSolexa_with_or_without_equalization(dirName,readLength,'O7')
createUniqueReads_12SamplesSolexa_with_or_without_equalization(dirName,readLength,'O10')
createUniqueReads_12SamplesSolexa_with_or_without_equalization(dirName,readLength,'S7')
createUniqueReads_12SamplesSolexa_with_or_without_equalization(dirName,readLength,'S10')
createUniqueReads_12SamplesSolexa_with_or_without_equalization(dirName,readLength,'M1')
createUniqueReads_12SamplesSolexa_with_or_without_equalization(dirName,readLength,'M2')
createUniqueReads_12SamplesSolexa_with_or_without_equalization(dirName,readLength,'M3')
createUniqueReads_12SamplesSolexa_with_or_without_equalization(dirName,readLength,'M4')
createUniqueReads_12SamplesSolexa_with_or_without_equalization(dirName,readLength,'S1')
createUniqueReads_12SamplesSolexa_with_or_without_equalization(dirName,readLength,'S2')
createUniqueReads_12SamplesSolexa_with_or_without_equalization(dirName,readLength,'S3')
createUniqueReads_12SamplesSolexa_with_or_without_equalization(dirName,readLength,'S4')

readLength = 100;
createUniqueReads_12SamplesSolexa_with_or_without_equalization(dirName,readLength,'O7')
createUniqueReads_12SamplesSolexa_with_or_without_equalization(dirName,readLength,'O10')
createUniqueReads_12SamplesSolexa_with_or_without_equalization(dirName,readLength,'S7')
createUniqueReads_12SamplesSolexa_with_or_without_equalization(dirName,readLength,'S10')
createUniqueReads_12SamplesSolexa_with_or_without_equalization(dirName,readLength,'M1')
createUniqueReads_12SamplesSolexa_with_or_without_equalization(dirName,readLength,'M2')
createUniqueReads_12SamplesSolexa_with_or_without_equalization(dirName,readLength,'M3')
createUniqueReads_12SamplesSolexa_with_or_without_equalization(dirName,readLength,'M4')
createUniqueReads_12SamplesSolexa_with_or_without_equalization(dirName,readLength,'S1')
createUniqueReads_12SamplesSolexa_with_or_without_equalization(dirName,readLength,'S2')
createUniqueReads_12SamplesSolexa_with_or_without_equalization(dirName,readLength,'S3')
createUniqueReads_12SamplesSolexa_with_or_without_equalization(dirName,readLength,'S4')

unix('mkdir /home/csfaculty/shental/CS/BAC/12Samples/Solexa/data/with_or_without_equalization/with_equalization;mkdir /home/csfaculty/shental/CS/BAC/12Samples/Solexa/data/with_or_without_equalization/without_equalization')

% move with and without to their directories


%%%%%%%%%%%%%%%%%%55
  
end

if 1==2 % prepare the run files equalization
  % without equalization
  testData12SamplesSolexa_with_or_without_equalization({'O7','O10','S7','S10','M1','M2','M3','M4','S1','S2','S3','S4'},'without_eq',50,0)
  testData12SamplesSolexa_with_or_without_equalization({'O7','O10','S7','S10','M1','M2','M3','M4','S1','S2','S3','S4'},'without_eq',100,0)
  
  % with equalization
  testData12SamplesSolexa_with_or_without_equalization({'O7','O10','S7','S10','M1','M2','M3','M4','S1','S2','S3','S4'},'with_eq',50,1)
  testData12SamplesSolexa_with_or_without_equalization({'O7','O10','S7','S10','M1','M2','M3','M4','S1','S2','S3','S4'},'with_eq',100,1)
  
end

if 1==2
  createFilesFor_testData12Samples_with_or_without_equalization(0,50)
  createFilesFor_testData12Samples_with_or_without_equalization(0,100)
  createFilesFor_testData12Samples_with_or_without_equalization(1,50)
  createFilesFor_testData12Samples_with_or_without_equalization(1,100)
end

if 1==2
   mkdir resSolexa_with_or_without_equalization
  % rm -f list_with_or_without_equalization;find sol_with* -name "*.mat" >> list_with_or_without_equalization
  % basj
  % for i in `cat list_with_or_without_equalization`; do cp $i resSolexa_with_or_without_equalization;done
end


%./run_no_correction_without_eq_50.txt;./run_with_correction_without_eq_50.txt;./run_no_correction_without_eq_100.txt;./run_with_correction_without_eq_100.txt;./run_no_correction_with_eq_50.txt;./run_with_correction_with_eq_50.txt;
%run_no_correction_with_eq_100.txt;
%run_with_correction_with_eq_100.txt;
