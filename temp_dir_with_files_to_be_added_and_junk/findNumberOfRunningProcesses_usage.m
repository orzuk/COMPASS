function findNumberOfRunningProcesses_usage(a)
%rmpath(genpath('/home/csfaculty/shental/Weizmann/ph2users/fenoam/'));rmpath('/home/csfaculty/shental/matlab');cd ~/CS/mFilesBAC;addpath('~/CS/mFilesBAC'); mcc -v -m findNumberOfRunningProcessesForCorrection.m -d /homes/csfaculty/shental/CS/BAC/cvx

flag = 1;
while flag==1
  [junk,val] = unix(['echo "start\n";echo "hour SUSP" ; bjobs -u orzuk | grep hour | grep SUSP | wc ; ',...
                     'echo "hour PEND" ; bjobs -u orzuk | grep hour | grep PEND | wc ; ',...
                     'echo "hour RUN" ; bjobs -u orzuk | grep hour | grep RUN | wc ;',...
                     ' echo "week SUSP" ; bjobs -u orzuk | grep week | grep SUSP | wc ; ',...
                     'echo "week PEND" ; bjobs -u orzuk | grep week | grep PEND | wc ; ',...,
                     'echo "week RUN" ; bjobs -u orzuk | grep week | grep RUN | wc ;']);

  a = findstr(val,'start');
  val = val(a(1)+7:end);
  a = find(val==10);
  
  k = 1;
  for i=1:length(a)-1
    b = str2num(val(a(i)+1:a(i+1)-1));
    if ~isempty(b) && isnumeric(b)
      v(k,:) = b;
      k = k+1;
    end
    %
  end
  
  if sum(v(4:6,1))>10 | sum(v(1:3,1))>500
    disp('pause till are done')
    pause(60)
  else
    flag = 0;
    disp('running again')
  end
  
  
end
