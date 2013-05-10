if 1==2 % creating reads

% create reads
clear

dirName = 'with_eq_unicov_mar21';
typeName = '_unicov_corrfreq';

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
  testData12Samples_with_eq_unicov_corrfreq({'S1','O7','O10','S7','S10','M1','M2','M3','M4','S2','S3','S4'},'with_eq_unicov_mar21','with_eq_unicov_mar21')
end

if 1==2
  % create files for run
  createFilesFor_testData12Samples_with_eq_unicov('with_eq_unicov_mar21')
end


if 1==2
  mkdir resSolexa_with_eq_unicov_mar21
  % rm -f list_with_eq_unicov_mar21;find sol_with_eq_unicov_mar21* -name "*.mat" >> list_with_eq_unicov_mar21
  % basj
  % for i in `cat list_with_eq_unicov_mar21`; do cp $i resSolexa_with_eq_unicov_mar21;done
end

%%%%%%%%%%%%%%%%555
% test
clear
dirName = 'with_eq_unicov_first_try';
typeName = '_unicov_corrfreq';

readLength = 100;
createUniqueReads_12SamplesSolexa_with_eq_unicov_corrfreq(dirName,readLength,'S1',typeName)

testData12Samples_with_eq_unicov_corrfreq({'S1'},'with_eq_unicov_first_try','with_eq_unicov_first_try')

createFilesFor_testData12Samples_with_eq_unicov('with_eq_unicov_first_try')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%run again with former db

%datasetName = 'primers750_primer_tail';
datasetName = 'primers750_primer_tail_before_right_edge_correction';
testData12Samples_with_eq_unicov_corrfreq_again({'S1','O7','O10','S7','S10','M1','M2','M3','M4','S2','S3','S4'},'test_unicov_mar22_with_primer_tail_before_right_edge_correction','with_eq_unicov_mar21',datasetName)
createFilesFor_testData12Samples_with_eq_unicov('test_unicov_mar22_with_primer_tail_before_right_edge_correction')

dirName = 'test_unicov_mar22_with_primer_tail_before_right_edge_correction';
unix(['cd /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/;mkdir resSolexa_',dirName,';','rm -f /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/list_',dirName,';find sol_',dirName,'*  -name "*.mat" > list_',dirName])
% bash
w = ['for i in `cat list_',dirName,'`; do cp $i resSolexa_',dirName,';done'];


datasetName = 'primers750_before_right_edge_correction';
testData12Samples_with_eq_unicov_corrfreq_again({'S1','O7','O10','S7','S10','M1','M2','M3','M4','S2','S3','S4'},'test_unicov_mar22_primers750_before_right_edge_correction','with_eq_unicov_mar21',datasetName)
createFilesFor_testData12Samples_with_eq_unicov('test_unicov_mar22_primers750_before_right_edge_correction')

dirName = 'test_unicov_mar22_primers750_before_right_edge_correction';
unix(['cd /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/;mkdir resSolexa_',dirName,';','rm -f /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/list_',dirName,';find sol_',dirName,'*  -name "*.mat" > list_',dirName])
% bash
w = ['for i in `cat list_',dirName,'`; do cp $i resSolexa_',dirName,';done'];

