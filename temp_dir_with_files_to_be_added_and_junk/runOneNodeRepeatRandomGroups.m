function runOneNodeRepeatRandomGroups(fileToLoad,l,runFileName)

fid = fopen([runFileName],'w');
while fid<0
  pause(5)
  fid = fopen([runFileName],'w');
end

fprintf(fid,'#!/bin/bash \n');
  
fprintf(fid,'TMPDIR=/tmp/${RANDOM}_${RANDOM}\n');
fprintf(fid,'mkdir $TMPDIR\n');
fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
fprintf(fid,'%s/CS/BAC/cvx/run_distributeBAC_repeatRandomGroups.sh /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b %s %s\n',...
        fileToLoad,num2str(l));
fprintf(fid,'\n');
fprintf(fid,'rm -rf $TMPDIR \n');
fclose(fid);


