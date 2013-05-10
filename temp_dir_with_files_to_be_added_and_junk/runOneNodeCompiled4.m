function runOneNodeCompiled4(allocateToOne,basicAllocationFileName,mFileName,saveFileName,userDir,marilynNumber,queueName)
%keyboard
% the same as runOneNodeCompiled3 but runs runOneNodeForCompilation4
disp('runOneNodeCompiled4.m  ')

if ~exist('queueName')
  queueName = 'priority';
end

if isempty(findstr(userDir,'/seq/'))
  if exist('marilynNumber')
    %unix(['qsub  -q all.q@marilyn-',num2str(marilynNumber),...
    %      ' -o ',mFileName,'.o ',...
    %      ' -e ',mFileName,'.e ',...      
    %      userDir,'/CS/BAC/cvx/run_runOneNodeForCompilation4_WIS.sh /storage/MATLAB/R2010b/ ',num2str(allocateToOne(1)),' ',num2str(allocateToOne(end)) ,' ',basicAllocationFileName,' ',saveFileName,' ',userDir]);
  
  else
    
    
  end
else % Br    
  fid = fopen([mFileName],'w');
  fprintf(fid,'#!/bin/bash \n');
  
  fprintf(fid,'TMPDIR=/tmp/${RANDOM}_${RANDOM}\n');
  fprintf(fid,'mkdir $TMPDIR\n');
  fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
  fprintf(fid,'%s/CS/BAC/cvx/run_runOneNodeForCompilation4_Br.sh /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b %s %s %s %s %s %s\n',userDir,num2str(allocateToOne(1)),' ',num2str(allocateToOne(end)) ,basicAllocationFileName,saveFileName,userDir);

  fprintf(fid,'rm -rf $TMPDIR \n');
  fclose(fid);
end
% requeues for specific exit codes
% -Q "1 134"

%fprintf('%s %s %s %s %s %s\n')


