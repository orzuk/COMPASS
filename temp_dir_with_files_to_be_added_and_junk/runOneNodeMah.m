function runOneNodeMah(allocateToOne,basicAllocationFileName,mFileName,saveFileName,userDir,marilynNumber,queueName)
%keyboard


if ~exist('queueName')
  queueName = 'priority';
end
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
fprintf(fid,'%s/CS/BAC/cvx/run_ResultsDistance_Br.sh /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b %s %s %s %s %s %s\n',...
        userDir, ...
        num2str(allocateToOne(1)),' ',num2str(allocateToOne(end)) ,'a',basicAllocationFileName,saveFileName);
% '1' is just replacemnet for auxData needed in solveForGroup3                                    
fprintf(fid,'\n');
fprintf(fid,'rm -rf $TMPDIR \n');
fclose(fid);


