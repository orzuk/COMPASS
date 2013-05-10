function [X,sumRelevantReads]=iterateParallelDistributedSequenceFiles3(uniqueReads,uniqueReads_length,curr_kp,basicSeqNameDir,basicSeqKey,parallelType,basicSaveName,userDir,runFileName,auxData)
%keyboard
% deleted the following input as they appear in auxData: readLength,numProcessors,,queueName

disp('assumes that auxData does not change in iterateParallelDistributedSequenceFiles3 since it is saved to the file basicAllocationFileName. iterateParallelDistributedSequenceFiles3.m')

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



part = 1:auxData.groupSize:n;
part(end) = n+1;

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
    w = ['matlabpool open local ',num2str(currNumProcessors),';parfor i=1:length(tmpInd),i,[xx,ss]=solveForGroup3(i,i,auxData,[],[],tmpInd,basicSeqNameDir,basicSeqKey,readLength,uniqueReads,uniqueReads_length);x{i} = xx{i};sumRelevantReads(i) = ss(i);;end;matlabpool close'];
    eval(w);
    
  %keyboard
  else % single
    %keyboard 
    for i=1:length(tmpInd)
      [res.x,res.sumRelevantReads]=solveForGroup3(i,i,auxData,[],[],tmpInd,basicSeqNameDir,basicSeqKey,readLength,uniqueReads,uniqueReads_length);
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
      runOneNodeBasedSolveForGroup3(allocateToOne,basicAllocationFileName,mFileName,saveFileName,userDir,[],auxData.queueName);
    end
    %keyboard
    % chmod
    unix(['chmod 0700 ',runFileName,'_partDist_*']);
    
    
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
    fprintf(fid,'bsub   -q  %s -J "%s[1-%d]" -o "%s_partDist_Z%s.o" ~/runInput_donotdelete %s_partDist_%s;\n',auxData.queueName,tmpRunFileName,length(partDist)-1,runFileName,'%I',runFileName,'\$LSB_JOBINDEX');
    fprintf(fid,'rm -rf $TMPDIR \n');
    fclose(fid)
    %keyboard
    unix(['chmod 0700 ',runFileName,'_shell.sh;'])
    
    
    JobList{1} = tmpRunFileName;
    
    % submit these and test if were submitted
    %keyboard
    wereTheySubmittedFlag = 0;
    while wereTheySubmittedFlag==0
      unix([runFileName,'_shell.sh']);
      % check whether these jobs were indeed submitted
      pause(60)
      wereTheySubmittedFlag = testSubmitted(JobList{1});
    end
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% end running in batch
    
    
    numSuccessfullyCompleted = zeros(length(partDist)-1,1);
    while length(find(numSuccessfullyCompleted))<length(partDist)-1
      
      % check whether it got thrown away - leaving .o files 
      logSuccessfulFileFlag = 0;
      while logSuccessfulFileFlag==0
        unix(['grep -l "Successfully completed" ',runFileName,'_partDist_Z*.o | cut -f2 -d"Z" | cut -f1 -d"." >',runFileName,'_logOfSuccessful.txt;sleep   2;']);
        [j1,j2] = unix(['ls ',runFileName,'_logOfSuccessful.txt']);
        if j1==0
          disp([runFileName,'_logOfSuccessful.txt was found']);
          logSuccessfulFileFlag = 1;
        end
        pause(2)
      end
      
      logUnsuccessful1 = 0;
      while logUnsuccessful1==0
        unix(['grep -l "exit code" ',runFileName,'_partDist_Z*.o | cut -f2 -d"Z" | cut -f1 -d"." >',runFileName,'_logOfUnsuccessful1.txt;sleep 2;']);
        [j1,j2] = unix(['ls ',runFileName,'_logOfUnsuccessful1.txt;']);
        if j1==0
          disp([runFileName,'_logOfUnsuccessful1.txt was found']);
          logUnsuccessful1 = 1;
        end
        pause(2)
        
      end
      
      logUnsuccessful2 = 0;
      while logUnsuccessful2==0
        unix(['grep -l "fatalexit" ',runFileName,'_partDist_Z*.o | cut -f2 -d"Z" | cut -f1 -d"." >',runFileName,'_logOfUnsuccessful2.txt;sleep 2;']);
        [j1,j2] = unix(['ls ',runFileName,'_logOfUnsuccessful2.txt;']);
        if j1==0
          disp([runFileName,'_logOfUnsuccessful2.txt was found']);
          logUnsuccessful2 = 1;
        end
        pause(2)
      end
      
      logUnsuccessful3 = 0;
      while logUnsuccessful3==0
        unix(['grep -l "TERM_RUNLIMIT" ',runFileName,'_partDist_Z*.o | cut -f2 -d"Z" | cut -f1 -d"." >',runFileName,'_logOfUnsuccessful3.txt;']);
        [j1,j2] = unix(['ls ',runFileName,'_logOfUnsuccessful3.txt']);
        if j1==0
          disp([runFileName,'_logOfUnsuccessful3.txt was found']);
          logUnsuccessful3 = 1;
        end
        pause(2)
      end

      % load these
      listOfsuccessful = load([runFileName,'_logOfSuccessful.txt']);
      listOfProblematicBasedExitCode1 = load([runFileName,'_logOfUnsuccessful1.txt']);
      listOfProblematicBasedExitCode2 = load([runFileName,'_logOfUnsuccessful2.txt']);
      listOfProblematicBasedExitCode3 = load([runFileName,'_logOfUnsuccessful3.txt']);
      listOfProblematicBasedExitCode = [listOfProblematicBasedExitCode1;listOfProblematicBasedExitCode2;listOfProblematicBasedExitCode3];
      
      % did not run
      listOfProblematicDidNotRunButAppearInSuccesfulList = [];
      for j=listOfsuccessful'
        dr = dir([basicAllocationFileName,'_partDist_',num2str(j),'_res.mat']);
        if length(dr)==1
          numSuccessfullyCompleted(j) = 1;
        else
          listOfProblematicDidNotRunButAppearInSuccesfulList = [listOfProblematicDidNotRunButAppearInSuccesfulList,j];
        end
      end
      
      
      listOfProblematic = [listOfProblematicBasedExitCode;listOfProblematicDidNotRunButAppearInSuccesfulList'];
      listOfProblematic = unique(listOfProblematic);
      
      if ~isempty(listOfProblematic)
        pause(2)
        disp('running problematic. iterateParallelDistributedSequenceFiles3.m')
        % delete the .o files and create names
        
        currJobList = runProblematic3(listOfProblematic,runFileName,tmpRunFileName,auxData);
        for jj=1:length(currJobList)
          JobList{end+1} = currJobList{jj};
        end
        
        % delete the *.txt files - so that it does not run again
        unix(['rm -f ',runFileName,'_logOf*.txt']);
      
      else % wait for results/problems
        pause(2)
      end
      
      %%%%%%%%%%% end % check whether it got thrown away - leaving .0 files 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
      
      
      
      
      
      % check the number of res files - if ok - break
      resNames = dir([basicAllocationFileName,'_partDist_*_res.mat']);
      if   length(resNames)==length(partDist)-1
        disp('did all parts. iterateParallelDistributedSequenceFiles3.m')
        break
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%
      
      
      
      
      
      % run processes that started - i.e. were given a job number but were forgotten on the side
      
      % find the number of active processes - if zero processes are running and no res file exist - run them again
      
      processesAreStillRunningFlag = 0;
      jj = 1;
      while processesAreStillRunningFlag==0 && jj<=length(JobList)
        [j1,j2] = unix(['bjobs -J  ',JobList{jj},' >',runFileName,'_',num2str(jj),'_checkIfFininished.txt',';grep "is not found" ',runFileName,'_',num2str(jj),'_checkIfFininished.txt',' > ',runFileName,'_',num2str(jj),'_checkIfFininished.txt']);
        if isempty(findstr(j2,'is not found'))
          processesAreStillRunningFlag = 1;
          disp('at least one processes if still on')
        end
        jj = jj+1;
      end
     
      if processesAreStillRunningFlag==0
        disp('all processes are off. checking if all res exist.')
        resNames = dir([basicAllocationFileName,'_partDist_*_res.mat']);
        if   length(resNames)==length(partDist)-1
          disp('did all parts. iterateParallelDistributedSequenceFiles3.m')
          break
        else
          disp('running processes that were forgotten along the way. iterateParallelDistributedSequenceFiles3.m')
          listOfProblematic = find(numSuccessfullyCompleted==0);
          JobList = [];
          JobList = runProblematic3(listOfProblematic,runFileName,tmpRunFileName,auxData);
          
        end
      end
      %%%%%%%%%%%% end  run processes that started - i.e. were given a job number but were forgotten on the side
     
      
      
      
      pause(60)
      
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
  unix(['rm -f ',auxData.currDir,'/core*']);
  unix(['rm -f ',runFileName,'_*'])
  unix(['rm -f ',basicSaveName,'_*'])
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
    disp('problem. iterateParallelDistributedSequenceFiles3.m')
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






  
