function findNumberOfRunningProcessesFourth(userDir,fileNameForFindNumber)

%rmpath(genpath('/home/csfaculty/shental/Weizmann/ph2users/fenoam/'));rmpath('/home/csfaculty/shental/matlab');cd ~/CS/mFilesBAC;addpath('~/CS/mFilesBAC'); mcc -v -m findNumberOfRunningProcessesFourth.m -d /homes/csfaculty/shental/CS/BAC/cvx

load([userDir,'/CS/BAC/cvx/',fileNameForFindNumber]); % hold a auxDataFileNumber



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
  
  num_hour_susp = v(1,1);
  num_hour_pend = v(2,1);
  num_hour_run = v(3,1);
  num_week_susp = v(4,1);
  num_week_pend = v(5,1);
  num_week_run = v(6,1);
  
  
  
  %if num_hour_susp+num_hour_pend+num_hour_run>maxNumberInHour | ...
  %   num_week_susp+num_week_pend+num_week_run>maxNumberInWeek | ...
  %   num_week_pend>maxNumWeekPending
  if   num_hour_susp+num_hour_pend>auxDataFileNumber.maxNumberHourPending |...
       num_week_pend > auxDataFileNumber.maxNumberWeekPending | ...
       num_hour_run > auxDataFileNumber.maxNumberHourRunning
    
    disp('pause till are done')
    pause(auxDataFileNumber.pauseTimeSeconds)
    
  else
    flag = 0;
    disp('running again')
  end
  
  
end
