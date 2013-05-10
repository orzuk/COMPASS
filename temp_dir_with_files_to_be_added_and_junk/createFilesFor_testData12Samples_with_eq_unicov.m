function createFilesFor_testData12Samples_with_eq_unicov(datasetName)

readLength = 100;

currFileName = ['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/run_no_correction_',datasetName,'.txt'];

fid = fopen(currFileName,'w');
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
fprintf(fid,'mkdir $TMPDIR\n');
fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
fprintf(fid,'\n');  

sampleName = {'S1','O7','O10','S7','S10','M1','M2','M3','M4','S2','S3','S4'};
for i=1:length(sampleName)
    fprintf(fid,'/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_%s_%s/run_sol_%s_%s_noCorrection_%d\n',datasetName,sampleName{i},datasetName,sampleName{i},readLength);
end

fprintf(fid,'rm -fr $TMPDIR\n');
fclose(fid);
unix(['chmod 0700 ',currFileName]);




