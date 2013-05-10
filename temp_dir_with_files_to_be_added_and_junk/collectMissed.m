function [fileNumbersMissed,fileNumbersDidNotRun]=collectMissed(basicPrefix,outFileDirName,iterationNumber)
%userDir = getuserdir;
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
keyboard
unix(['ls -l ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',basicPrefix,'_bac_dist_flag_*.e | tr -s " " |cut -f5,9  -d" ">  ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','problamatic_',basicPrefix,'_iter_',num2str(iterationNumber)]);

fid = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','problamatic_',basicPrefix,'_iter_',num2str(iterationNumber)],'r');
list = textscan(fid,'%d %s');
fclose(fid);

fclose('all');
k = 1;
fid = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','listOfR_rerun_',num2str(iterationNumber),'_',basicPrefix,'_',num2str(k)],'w');
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
fprintf(fid,'mkdir $TMPDIR\n');
fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
fprintf(fid,'\n');

fileNumbersMissed = k;
for i=1:size(list{1,1},1)
  if list{1,1}(i)>0
    
    if mod(k,10)==1
      if exist('fid')
        fprintf(fid,'rm -fr $TMPDIR\n');
        fclose(fid);
        clear fid
      end
      
      fileNumbersMissed = [fileNumbersMissed,k];
      fid = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','listOfR_rerun_',num2str(iterationNumber),'_',basicPrefix,'_',num2str(k)],'w');
      fprintf(fid,'#!/bin/bash\n');
      fprintf(fid,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
      fprintf(fid,'mkdir $TMPDIR\n');
      fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
      fprintf(fid,'\n');
      
    end
    
    name = list{1,2}{i};
    
    % del .e *.o
    name(end-1:end) = [];
    fprintf(fid,'rm -f %s %s\n',[name,'.e'],[name,'.o']); % del .e .o 
    
    
    
    %del their position
    a = find(name=='/');
    pre = name(1:a(end-2));
    post = name(a(end)+1:end);
    fprintf(fid,'rm -f %s/* \n',[pre,'tmpRuns/',outFileDirName,'/',post])
    
     
    % rerun
    a = find(name=='/');
    post = name(a(end)+1:end);
    pre = name(1:a(end));
    fprintf(fid,'%s\n',[pre,'run_',post]);
 
    fprintf(fid,'\n');
    
    k = k+1;
  end % list{1}
end
if exist('fid')
  fprintf(fid,'rm -fr $TMPDIR\n');
  fclose(fid);
  clear fid
end


% del *.e *.o



unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','listOfR_rerun_',num2str(iterationNumber),'_',basicPrefix,'_*']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect those that did not run
load([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',basicPrefix,'_list'])


didNotRunList = [];
for i=1:size(totalListOfNames,2)
  flagFound = 0;
  for j=1:size(list{1,2},1)
    nameList = list{1,2}{j};
    nameList(end-1:end) = [];
    a = find(nameList=='/');
    nameList = nameList(a(end)+1:end);
    if strcmp(totalListOfNames{i},nameList);
      flagFound = 1;
    end
    
  end
  
  if flagFound==0
    didNotRunList{end+1} = totalListOfNames{i};
  end
  
end

% run those that did not run at all

if length(didNotRunList)==0
  fileNumbersDidNotRun = 0;
  return
end


name = list{1,2}{1};
a = find(name=='/');
pre = name(1:a(end));

fclose('all');
k = 1;
fid = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','listOfR_didNotRun_',num2str(iterationNumber),'_',basicPrefix,'_',num2str(k)],'w');
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
fprintf(fid,'mkdir $TMPDIR\n');
fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
fprintf(fid,'\n');

fileNumbersDidNotRun = k;
for i=1:length(didNotRunList)
  
  if mod(k,10)==1
      if exist('fid')
        fprintf(fid,'rm -fr $TMPDIR\n');
        fclose(fid);
        clear fid
      end
      
      fileNumbersDidNotRun = [fileNumbersDidNotRun,k];
      fid = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','listOfR_didNotRun_',num2str(iterationNumber),'_',basicPrefix,'_',num2str(k)],'w');
      fprintf(fid,'#!/bin/bash\n');
      fprintf(fid,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
      fprintf(fid,'mkdir $TMPDIR\n');
      fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
      fprintf(fid,'\n');
      
  end
  
  % run
 
  fprintf(fid,'%s\n',[pre,'run_',didNotRunList{i}]);
  fprintf(fid,'\n');
    
  k = k+1;
  
end
if exist('fid')
  fprintf(fid,'rm -fr $TMPDIR\n');
  fclose(fid);
  clear fid
end


unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','listOfR_didNotRun_',num2str(iterationNumber),'_',basicPrefix,'_*']);
