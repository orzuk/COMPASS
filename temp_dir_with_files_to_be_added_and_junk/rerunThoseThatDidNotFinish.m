function [fileNumbersMissed,fileNumbersDidNotRun]=rerunThoseThatDidNotFinish(basicPrefix,outFileDirName,iterationNumber)
%userDir = getuserdir;
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
keyboard
for i=1:size(list{1,1},1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect those that did not run
load([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',basicPrefix,'_list'])

unix(['ls -1 ',userDir,'/CS/BAC/',outFileDirName,'/*.mat > ',userDir,'/CS/BAC/',outFileDirName,'/done_',num2str(iterationNumber),'.txt']);
fid = fopen([userDir,'/CS/BAC/',outFileDirName,'/done_',num2str(iterationNumber),'.txt'],'r');
doneList = textscan(fid,'%s');
fclose(fid);


didNotRunList = [];
for i=1:size(totalListOfNames,2)
  flagFound = 0;
  for j=1:length(doneList{1})
    nameList = doneList{1}{j};
    nameList(end-3:end) = [];
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


pre = [userDir,'/CS/BAC/dataForSim/',outFileDirName,'/'];

fclose('all');
k = 1;
fid = fopen([pre,'listOfR_didNotRun_',num2str(iterationNumber),'_',basicPrefix,'_',num2str(k)],'w');
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
