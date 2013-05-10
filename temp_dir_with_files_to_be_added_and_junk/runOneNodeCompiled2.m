function runOneNodeCompiled2(allocateToOne,basicAllocationFileName,mFileName,saveFileName,userDir,marilynNumber,queueName)
%keyboard

disp('runOneNodeCompiled2. The only difference is that it runs: run_runOneNodeForCompilation_WIS2.sh instead of run_runOneNodeForCompilation_WIS.sh')

if ~exist('queueName')
  queueName = 'priority';
end

if isempty(findstr(userDir,'/seq/'))
  if exist('marilynNumber')
    unix(['qsub  -q all.q@marilyn-',num2str(marilynNumber),...
          ' -o ',mFileName,'.o ',...
          ' -e ',mFileName,'.e ',...      
          userDir,'/CS/BAC/cvx/run_runOneNodeForCompilation2_WIS.sh /storage/MATLAB/R2010b/ ',num2str(allocateToOne(1)),' ',num2str(allocateToOne(end)) ,' ',basicAllocationFileName,' ',saveFileName,' ',userDir]);
  
  else
    
    
  end
else % Br
  unix(['bsub   -q ',queueName,...
          ' -M 33554432 ',...
          ' -o ',mFileName,'.o ',...
          ' -e ',mFileName,'.e ',...      
          userDir,'/CS/BAC/cvx/run_runOneNodeForCompilation2_Br.sh /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b ',num2str(allocateToOne(1)),' ',num2str(allocateToOne(end)) ,' ',basicAllocationFileName,' ',saveFileName,' ',userDir]);
end
% requeues for specific exit codes
% -Q "1 134"