function runOneNodeBasedSolveForGroupOr(allocateToOne,basicAllocationFileName,mFileName,saveFileName,userDir,marilynNumber,queueName)
%keyboard


if ~exist('queueName')
  queueName = 'priority';
end

if isempty(findstr(userDir,'/seq/'))
  if exist('marilynNumber')
    %unix(['qsub  -q all.q@marilyn-',num2str(marilynNumber),...
    %      ' -o ',mFileName,'.o ',...
    %      ' -e ',mFileName,'.e ',...      
    %      userDir,'/CS/BAC/cvx/run_runOneNodeForCompilation2_WIS.sh /storage/MATLAB/R2010b/ ',num2str(allocateToOne(1)),' ',num2str(allocateToOne(end)) ,' ',basicAllocationFileName,' ',saveFileName,' ',userDir]);
  
  else
    
    
  end
else % Br    
  fid = fopen([mFileName],'w');
  while fid<0
    pause(5)
    fid = fopen([mFileName],'w');
  end
  fprintf(fid,'#!/bin/bash \n');
  
  fprintf(fid,'TMPDIR=/tmp/${RANDOM}_${RANDOM}\n');
  fprintf(fid,'mkdir $TMPDIR\n');
  fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
  fprintf(fid,'%s/CS/BAC/cvx/run_solveForGroupOr_Br.sh /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b %s %s %s %s %s %s\n',...
                    userDir, ...
num2str(allocateToOne(1)),' ',num2str(allocateToOne(end)) ,'a',basicAllocationFileName,saveFileName);
                         % '1' is just replacemnet for auxData needed in solveForGroup3                                    
  fprintf(fid,'\n');
  fprintf(fid,'rm -rf $TMPDIR \n');
  fclose(fid);
end
% requeues for specific exit codes
% -Q "1 134"

%fprintf('%s %s %s %s %s %s\n')

