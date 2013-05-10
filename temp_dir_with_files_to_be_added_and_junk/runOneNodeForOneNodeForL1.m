function runOneNodeForOneNodeForL1(basicAllocationFileName,mFileName,userDir)
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
fprintf(fid,'%s/CS/BAC/cvx/run_oneNodeForL1OrFourth.sh /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b %s\n',userDir,basicAllocationFileName);
fprintf(fid,'\n');
fprintf(fid,'rm -rf $TMPDIR \n');
fclose(fid);


