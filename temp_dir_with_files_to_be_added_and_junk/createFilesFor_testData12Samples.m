function createFilesFor_testData12Samples


% no correction
currFileName = ['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/run_no_correction_sliding50.txt'];
fid = fopen(currFileName,'w');
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
fprintf(fid,'mkdir $TMPDIR\n');
fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
fprintf(fid,'\n');  

sampleName = {'O7','O10','S7','S10','M1','M2','M3','M4','S1','S2','S3','S4'};
for i=1:length(sampleName)
  fprintf(fid,'/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_WITH_REVERSE_sliding50_%s/run_sol_WITH_REVERSE_sliding50_%s_noCorrection_50\n',sampleName{i},sampleName{i});
end
fprintf(fid,'rm -fr $TMPDIR\n');
fclose(fid);
unix(['chmod 0700 ',currFileName]);



% with correction
currFileName = ['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/run_with_correction_sliding50.txt'];
fid = fopen(currFileName,'w');
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
fprintf(fid,'mkdir $TMPDIR\n');
fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
fprintf(fid,'\n');  

sampleName = {'O7','O10','S7','S10','M1','M2','M3','M4','S1','S2','S3','S4'};
for i=1:length(sampleName)
  fprintf(fid,'/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_WITH_REVERSE_sliding50_%s/run_sol_WITH_REVERSE_sliding50_%s_withCorrection_50\n',sampleName{i},sampleName{i});
end

fprintf(fid,'rm -fr $TMPDIR\n');
fclose(fid);
unix(['chmod 0700 ',currFileName]);

