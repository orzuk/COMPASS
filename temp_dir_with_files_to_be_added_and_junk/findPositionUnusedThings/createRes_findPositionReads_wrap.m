function createRes_findPositionReads_wrap(sampleName)


userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';

unix(['mkdir /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/findAlignment/tmpRuns/alignTmp_',sampleName])
unix(['mkdir /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/findAlignment/dataForSim/alignTmp_',sampleName])

currDir = [userDir,'/CS/BAC/findAlignment/dataForSim/alignTmp_',sampleName];
currFileName = [currDir,'/shell_',sampleName,'.sh'];

curr_fid = fopen(currFileName,'w');
fprintf(curr_fid,'#!/bin/bash\n');
fprintf(curr_fid,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
fprintf(curr_fid,'mkdir $TMPDIR\n');
fprintf(curr_fid,'export MCR_CACHE_ROOT=$TMPDIR\n');

master_queueName = 'week';
w1 = [sprintf(['cd ',currDir,';'])];
fprintf(curr_fid,'%s\n',w1);
w1 = ['bsub -c 1440 -q ',master_queueName,...
      ' -o ',userDirForFile,'/CS/BAC/dataForSim/',outFileDirName,'/',basicName,'.o ',...
      ' -e ',userDirForFile,'/CS/BAC/dataForSim/',outFileDirName,'/',basicName,'.e ',...
      '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/cvx/run_distributeBAC_generalOrFourth.sh  /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b ',auxDataFileName];
% deleted ' -M 33554432',
                  
fprintf(curr_fid,'%s\n',w1);
fprintf(curr_fid,'rm -fr $TMPDIR\n');
fclose(curr_fid);
clear curr_fid
                  



