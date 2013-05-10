function createFilesFor_testData12Samples_with_or_without_equalization(equalizationFlag,readLength)


% no correction
if equalizationFlag
  add = 'with_eq';
else
  add = 'without_eq';
end

currFileName = ['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/run_no_correction_',add,'_',num2str(readLength),'.txt'];

fid = fopen(currFileName,'w');
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
fprintf(fid,'mkdir $TMPDIR\n');
fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
fprintf(fid,'\n');  

sampleName = {'O7','O10','S7','S10','M1','M2','M3','M4','S1','S2','S3','S4'};
if equalizationFlag
  for i=1:length(sampleName)
    fprintf(fid,'/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_with_eq_%s/run_sol_with_eq_%s_noCorrection_%d\n',sampleName{i},sampleName{i},readLength);
  end
  
else
  for i=1:length(sampleName)
    fprintf(fid,'/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_without_eq_%s/run_sol_without_eq_%s_noCorrection_%d\n',sampleName{i},sampleName{i},readLength);
  end
end

fprintf(fid,'rm -fr $TMPDIR\n');
fclose(fid);
unix(['chmod 0700 ',currFileName]);



% with correction
currFileName = ['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/run_with_correction_',add,'_',num2str(readLength),'.txt'];

fid = fopen(currFileName,'w');
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
fprintf(fid,'mkdir $TMPDIR\n');
fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
fprintf(fid,'\n');  

sampleName = {'O7','O10','S7','S10','M1','M2','M3','M4','S1','S2','S3','S4'};
if equalizationFlag
  for i=1:length(sampleName)
    fprintf(fid,'/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_with_eq_%s/run_sol_with_eq_%s_withCorrection_%d\n',sampleName{i},sampleName{i},readLength);
  end
  
else
  for i=1:length(sampleName)
    fprintf(fid,'/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_without_eq_%s/run_sol_without_eq_%s_withCorrection_%d\n',sampleName{i},sampleName{i},readLength);
  end
end

fprintf(fid,'rm -fr $TMPDIR\n');
fclose(fid);
unix(['chmod 0700 ',currFileName]);

