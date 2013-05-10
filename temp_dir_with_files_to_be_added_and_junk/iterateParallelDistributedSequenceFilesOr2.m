function [X,sumRelevantReads]=iterateParallelDistributedSequenceFilesOr2(uniqueReads,uniqueReads_length,curr_kp,basicSeqNameDir,basicSeqKey,parallelType,basicSaveName,userDir,runFileName,auxData)
%keyboard
% deleted the following input as they appear in auxData: readLength,numProcessors,,queueName

disp('assumes that auxData does not change in iterateParallelDistributedSequenceFilesOr since it is saved to the file basicAllocationFileName. iterateParallelDistributedSequenceFilesOr.m')

readLength = auxData.readLength;
currNumProcessors = auxData.numProcessors;

n = length(curr_kp);

if auxData.keepOriginalOrderFlag==0
  currInd = randperm(n);
  disp('using permuted order')
else
  currInd = 1:n;
  disp('using the original order')
end

if auxData.forcePartFromOutsideFlag==0
  part = 1:auxData.groupSize:n;
  part(end) = n+1;
else
  disp('using order from outside. probably for repeatRandomGroups. iterateParallelDistributedSequenceFilesOr.m')
  part = auxData.partFromOutside;
end

%keyboard

tmpInd = cell(length(part)-1,1);
for i=1:length(part)-1
  tmpInd{i} = part(i):part(i+1)-1;
  tmpInd{i} = curr_kp(currInd(tmpInd{i}));
end



%keyboard
firstFlag = 1;
if strcmp(parallelType,'local')
  %keyboard
  %%%%%%%%
  % do not touch this
  if currNumProcessors>1 
    %keyboard
    
    w = ['matlabpool open local ',num2str(currNumProcessors),';parfor i=1:length(tmpInd),i,[xx,ss]=solveForGroupOr(i,i,auxData,[],[],tmpInd,basicSeqNameDir,basicSeqKey,readLength,uniqueReads,uniqueReads_length);x{i} = xx{i};sumRelevantReads(i) = ss(i);;end;matlabpool close'];
    eval(w);
    
  %keyboard
  else % single
    %keyboard 
    
    for i=1:length(tmpInd)
      [res.x,res.sumRelevantReads]=solveForGroupOr(i,i,auxData,[],[],tmpInd,basicSeqNameDir,basicSeqKey,readLength,uniqueReads,uniqueReads_length);
      x{i} = res.x{i};
      sumRelevantReads(i) = res.sumRelevantReads(i);
    end
  end % local
  
elseif strcmp(parallelType,'cluster')
  
  numParts = length(tmpInd);
  if numParts<=currNumProcessors
    currNumProcessors = numParts;
  end
  
  partDist = 1:round(numParts/currNumProcessors):numParts;
  partDist(end) = numParts+1;
  if length(partDist)==1
    partDist = [1,2];
  end
  
  
  % create the reads locations and sumRelevantReads
  basicAllocationFileName = [basicSaveName,'_basicAllocationData'];
  save(basicAllocationFileName,'readLength','basicSeqNameDir','basicSeqKey','uniqueReads','uniqueReads_length','tmpInd','auxData')
  
  if isempty(findstr(userDir,'/seq/')) % machon
    disp('uses Marilyn')
    cpuPerMachine = ceil(currNumProcessors/2);
    for j=1:length(partDist)-1
      if j<=cpuPerMachine
        marilynNumber = 1;
      elseif j>cpuPerMachine & j<=cpuPerMachine*2
        marilynNumber = 2;
      %elseif j>cpuPerMachine*2
      %  marilynNumber = 2;
      end
      
      allocateToOne = partDist(j):partDist(j+1)-1;
      mFileName = [basicAllocationFileName,'_partDist_',num2str(j)];
      saveFileName = [basicAllocationFileName,'_partDist_',num2str(j),'_res'];
      %runOneNode(allocateToOne,basicAllocationFileName,mFileName,saveFileName,userDir,marilynNumber);
      %keyboard
      runOneNodeCompiled2(allocateToOne,basicAllocationFileName,mFileName,saveFileName,userDir,marilynNumber,[]);
    end
    
    
    
    
  else % Br %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%keyboard
    % new way of submitting jobs - create the runFiles
    JobList = [];
    
    for j=1:length(partDist)-1
      allocateToOne = partDist(j):partDist(j+1)-1;
      mFileName = [runFileName,'_partDist_',num2str(j)];
      saveFileName = [basicAllocationFileName,'_partDist_',num2str(j),'_res'];
      runOneNodeBasedSolveForGroupOr(allocateToOne,basicAllocationFileName,mFileName,saveFileName,userDir,[],auxData.queueName);
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
      fprintf(fid,'bsub -W 1:30   -q  %s -W 1:59 -J "%s[1-%d]" -o "%s_partDist_Z%s.o" ~/runInput_donotdelete %s_partDist_%s;\n',auxData.queueName,tmpRunFileName,length(partDist)-1,runFileName,'%I',runFileName,'\$LSB_JOBINDEX');
      % -M 1000000
    else
      fprintf(fid,'bsub -q  %s -J "%s[1-%d]" -o "%s_partDist_Z%s.o" ~/runInput_donotdelete %s_partDist_%s;\n',auxData.queueName,tmpRunFileName,length(partDist)-1,runFileName,'%I',runFileName,'\$LSB_JOBINDEX');
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
    
    
    numSuccessfullyCompleted = zeros(length(partDist)-1,1);
    while length(find(numSuccessfullyCompleted))<length(partDist)-1
      
      
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
        dr = dir([basicAllocationFileName,'_partDist_',num2str(curr_i),'_res.mat']);
        if length(dr)==1
          numSuccessfullyCompleted(curr_i) = 1;
        else
          listOfProblematicDidNotRunButAppearInSuccesfulList = [listOfProblematicDidNotRunButAppearInSuccesfulList,curr_i];
        end
      end
      
      
      
      listOfProblematic = [listOfProblematicBasedExitCode;listOfProblematicDidNotRunButAppearInSuccesfulList'];
      listOfProblematic = unique(listOfProblematic);
      
       
      
      %[listOfProblematic1,numSuccessfullyCompleted1]=tmpTest(basicAllocationFileName,runFileName,partDist);
      
      %if ~isempty(listOfProblematic1)
      %  find(listOfProblematic1-listOfProblematic)
      %end
      %find(numSuccessfullyCompleted1-numSuccessfullyCompleted)
      %keyboard
      
      if ~isempty(listOfProblematic)
        %pause(2)
        disp('running problematic. iterateParallelDistributedSequenceFilesOr.m')
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
      resNames = dir([basicAllocationFileName,'_partDist_*_res.mat']);
      if   length(resNames)==length(partDist)-1
        disp('did all parts. iterateParallelDistributedSequenceFilesOr.m')
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
        resNames = dir([basicAllocationFileName,'_partDist_*_res.mat']);
        if   length(resNames)==length(partDist)-1
          disp('did all parts. iterateParallelDistributedSequenceFilesOr.m')
          break
        else
          disp('running processes that were forgotten along the way. iterateParallelDistributedSequenceFilesOr2.m')
          %keyboard
          listOfProblematic = find(numSuccessfullyCompleted==0);
          JobList = [];
          JobList = runProblematic3(listOfProblematic,runFileName,tmpRunFileName,auxData);
          
        end
      end
      %%%%%%%%%%%% end  run processes that started - i.e. were given a job number but were forgotten on the side
     
      
      
      
      pause(2)
      
    end %length(find(numSuccessfullyCompleted))<length(partDist)-1
    
    
  end %  Br 
  
  % check that all partDist have finished
  resNames = dir([basicAllocationFileName,'_partDist_*_res.mat']);
  while length(resNames)==0 || length(resNames)~=length(partDist)-1
    disp('checking whether it finished')
    pause(10)
    resNames = dir([basicAllocationFileName,'_partDist_*_res.mat']);
  end
  
  
  
  for j=1:length(partDist)-1
    allocateToOne = partDist(j):partDist(j+1)-1;
    res{j} = load([basicAllocationFileName,'_partDist_',num2str(j),'_res']);
    for i=allocateToOne
      x{i} = res{j}.x{i};
      sumRelevantReads(i) = res{j}.sumRelevantReads(i);
    end
  end
	
  % rm all irrelevant files
 
  %unix(['rm -f ',auxData.currDir,'/run_',auxData.tmpFileName,'_*']);
  
  %unix(['rm -f ',auxData.currDir,'/core*']);
  %unix(['rm -f ',runFileName,'_*'])
  %unix(['rm -f ',basicSaveName,'_*'])
  
  unix(['rm -f ',auxData.currDir,'/core*; ','rm -f ',runFileName,'_*; ','rm -f ',basicSaveName,'_*; ']);
   
  save([basicSaveName,'_tmpRes'],'x','sumRelevantReads')	
	    
end % strcmp(parallelType,'cluster')

% the results

size(x{1})


if size(x{1},2)==1 % standard
  X = zeros(1,n);
  for i=1:length(part)-1
    tmpInd = part(i):part(i+1)-1;
    tmpInd = currInd(tmpInd);
    if ~isempty(x{i})
      X(tmpInd) = x{i};
    end
  end
else % x_min x_max
  %keyboard
  if length(part)~=2 % should be a single one
    disp('problem. iterateParallelDistributedSequenceFilesOr.m')
    keyobard
  end
  %keyboard
  X = zeros(size(x{i}));
  for i=1:length(part)-1
    tmpInd = part(i):part(i+1)-1;
    tmpInd = currInd(tmpInd);
    if ~isempty(x{i})
      X(tmpInd,:) = x{i};
    end
  end
  
end




