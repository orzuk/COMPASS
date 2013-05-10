
% create reads
clear

dirName = 'opt4_th0_mar25_2012';
typeName = '_unifreq_corr_opt4_th0';
datasetName = 'primers750_primer_tail';

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

% prepare the run files equalization

testData12Samples_with_eq_unicov_corrfreq_again({'S1','O7','O10','S7','S10','M1','M2','M3','M4','S2','S3','S4'},['test_',dirName],dirName,datasetName)
createFilesFor_testData12Samples_with_eq_unicov(['test_',dirName])



unix(['cd /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/;mkdir resSolexa_test_',dirName,';','rm -f /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/list_test_',dirName,';find sol_test_',dirName,'*  -name "*.mat" > list_',dirName])
% bash
w = ['for i in `cat list_',dirName,'`; do cp $i resSolexa_test_',dirName,';done'];


