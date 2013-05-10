function doL1AsPartOf_distribute(basicAllocationFileName,runFileName,fileNameL1Res,auxData)
%keyboard
  
  % Br %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %keyboard
  % new way of submitting jobs - create the runFiles
  JobList = [];
  
  for j=1
    mFileName = [runFileName,'_partDist_',num2str(j)];
    runOneNodeForOneNodeForL1(basicAllocationFileName,mFileName,auxData.userDir);
  end
  %keyboard
  % chmod
  
  
  % run them in batch
  tmpRunFileName = runFileName;
  delSlash = find(tmpRunFileName=='/');
  tmpRunFileName = tmpRunFileName(delSlash(end)+1:end);
  %keyboard
  
  fid = fopen([runFileName,'_shell.sh'],'w');
  while fid<0 % error
    
    pause(5)
    fid = fopen([runFileName,'_shell.sh'],'w');
  end
  fprintf(fid,'#!/bin/bash\n');
  fprintf(fid,'TMPDIR=/tmp/${RANDOM}_${RANDOM}\n');
  fprintf(fid,'mkdir $TMPDIR\n');
  fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
  if strcmp(auxData.queueName,'hour')
    fprintf(fid,'bsub -W 1:30   -q  %s -W 1:59 -J "%s[1-%d]" -o "%s_partDist_Z%s.o" ~/runInput_donotdelete %s_partDist_%s;\n',auxData.queueName,tmpRunFileName,1,runFileName,'%I',runFileName,'\$LSB_JOBINDEX');
    % -M 1000000
  else
    fprintf(fid,'bsub -q  %s -J "%s[1-%d]" -o "%s_partDist_Z%s.o" ~/runInput_donotdelete %s_partDist_%s;\n',auxData.queueName,tmpRunFileName,1,runFileName,'%I',runFileName,'\$LSB_JOBINDEX');
    % -M 1000000 
  end
  
  fprintf(fid,'rm -rf $TMPDIR \n');
  fclose(fid)
  %keyboard
  unix(['chmod 0700 ',runFileName,'_partDist_*;','chmod 0700 ',runFileName,'_shell.sh;'])
  
  
  JobList{1} = tmpRunFileName;
  
  % submit these and test if were submitted
  %keyboard
  wereTheySubmittedFlag = 0;
  while wereTheySubmittedFlag==0
    unix([runFileName,'_shell.sh']);
    % check whether these jobs were indeed submitted
    pause(5)
    wereTheySubmittedFlag = testSubmitted(JobList{1});
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%% end running in batch
  
  
  numSuccessfullyCompleted = zeros(1,1);
  while length(find(numSuccessfullyCompleted))<1
    
    
    % check whether it got thrown away - leaving .o files 
    %keyboard
    
    
    
    %% new version  
    listOfFilesExitCode = grep('-l -s',{'exit code','fatalexit','TERM_RUNLIMIT'},[runFileName,'_partDist_Z*.o']);
    
    listOfProblematicBasedExitCode = zeros(1,length(listOfFilesExitCode));
    for i_lof=1:length(listOfFilesExitCode)
      index_lof_Z = findstr(listOfFilesExitCode{i_lof},'_Z');
      index_lof_o = findstr(listOfFilesExitCode{i_lof},'.o');
      listOfProblematicBasedExitCode(i_lof) = str2num(listOfFilesExitCode{i_lof}(index_lof_Z+2:index_lof_o-1));
    end
    listOfProblematicBasedExitCode = unique(listOfProblematicBasedExitCode);
    
    % did not run
    listOfFilesSuccessful = grep('-l -s',{'Successfully completed'},[runFileName,'_partDist_Z*.o']);
    listOfProblematicDidNotRunButAppearInSuccesfulList = [];
    for i_los=1:length(listOfFilesSuccessful)
      index_los_Z = findstr(listOfFilesSuccessful{i_los},'_Z');
      index_los_o = findstr(listOfFilesSuccessful{i_los},'.o');
      curr_i = str2num(listOfFilesSuccessful{i_los}(index_los_Z+2:index_los_o-1));
      dr = dir([fileNameL1Res,'.mat']);
      if length(dr)==1
        numSuccessfullyCompleted(curr_i) = 1;
      else
        listOfProblematicDidNotRunButAppearInSuccesfulList = [listOfProblematicDidNotRunButAppearInSuccesfulList,curr_i];
      end
    end
    
    
    
    listOfProblematic = [listOfProblematicBasedExitCode;listOfProblematicDidNotRunButAppearInSuccesfulList'];
    listOfProblematic = unique(listOfProblematic);
    
     %keyboard
    
    if ~isempty(listOfProblematic)
      %pause(2)
      disp('running problematic. doL1AsPartOf_distribute.m')
      % delete the .o files and create names
      
      currJobList = runProblematic3(listOfProblematic,runFileName,tmpRunFileName,auxData);
      for jj=1:length(currJobList)
        JobList{end+1} = currJobList{jj};
      end
      
      % delete the *.txt files - so that it does not run again
      %unix(['rm -f ',runFileName,'_logOf*.txt']);
      
    else % wait for results/problems
      pause(2)
    end
    
    %%%%%%%%%%% end % check whether it got thrown away - leaving .0 files 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    
    
    
    
    % check the number of res files - if ok - break
    resNames = dir([fileNameL1Res,'.mat']);
    if   length(resNames)==1
      disp('did all parts. doL1AsPartOf_distribute.m')
      break
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    % run processes that started - i.e. were given a job number but were forgotten on the side
    
    % find the number of active processes - if zero processes are running and no res file exist - run them again
        
    %keyboard
    
    processesAreStillRunningFlag = 0;
    wJobList = [];
    for jj=1:length(JobList)
      wJobList = [wJobList,'bjobs -J  ',JobList{jj},'| grep "is not found";'];
    end
    [j1,j2] = unix([wJobList]);
    if length(findstr(j2,'is not found'))<length(JobList)
      processesAreStillRunningFlag = 1;
    end
    
    
    if processesAreStillRunningFlag==0
      disp('all processes are off. checking if all res exist.')
      resNames = dir([fileNameL1Res,'.mat']);
      if   length(resNames)==1
        disp('did all parts. doL1AsPartOf_distrubte.m')
        break
      else
        disp('running processes that were forgotten along the way. doL1AsPartOf_distribute.m')
        %keyboard
        listOfProblematic = find(numSuccessfullyCompleted==0);
        JobList = [];
        JobList = runProblematic3(listOfProblematic,runFileName,tmpRunFileName,auxData);
        
      end
    end
    %%%%%%%%%%%% end  run processes that started - i.e. were given a job number but were forgotten on the side
    
    
    
    
    pause(2)
    
  end %length(find(numSuccessfullyCompleted))<length(partDist)-1
  
  
  %  Br 
  
  % check that all partDist have finished
  resNames = dir([fileNameL1Res,'.mat']);
  while length(resNames)==0 || length(resNames)~=1
    disp('checking whether it finished')
    pause(10)
    resNames = dir([fileNameL1Res,'.mat']);
  end
  
  
  unix(['rm -f ',auxData.currDir,'/core*; ','rm -f ',runFileName,'_*; ']);
   
	    
end % strcmp(parallelType,'cluster')


