function runOneNodeAlignment(mFileName,saveFileName,seq,indexFirst,indexLast,auxData)
%keyboard


% Br  
fid = fopen([mFileName],'w');
while fid<0
  pause(5)
  fid = fopen([mFileName],'w');
end
fprintf(fid,'#!/bin/bash \n');  
fprintf(fid,'TMPDIR=/tmp/${RANDOM}_${RANDOM}\n');
fprintf(fid,'mkdir $TMPDIR\n');
fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
fprintf(fid,'/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/cvx/run_myFindSeqInDB_generic.sh /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b %s %s %s %s %s %s %s %s %s\n',...
        auxData.basicName,auxData.basicDir,auxData.databaseName,auxData.blastPath,num2str(indexFirst),num2str(indexLast),saveFileName,auxData.dataDir,auxData.dataName);
% '1' is just replacemnet for auxData needed in solveForGroup3                                    
fprintf(fid,'\n');
fprintf(fid,'rm -rf $TMPDIR \n');
fclose(fid);


