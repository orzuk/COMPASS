function distribute_Mah(auxData)

%if isdeployed  %load a file which has a variable by the name auxData 
%  load(auxDataFileName)
%else
%  load(auxDataFileName)
%end

if ~isfield(auxData,'brFlag')
  auxData.brFlag = 0;
end

if auxData.brFlag
  userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
  disp('run in Br')
  parallelType = 'cluster';
else
  userDir = getuserdir;
  [junk systemName] = system('hostname');
  systemName = deblank(systemName);
  if ~isempty(findstr(userDir,'fenoam')) & ~strcmp(systemName,'complex1') %WIS
    parallelType = 'cluster';
    disp('runs in WIS')
  elseif ~isempty(findstr(userDir,'fenoam')) & strcmp(systemName,'complex1') %WIS
    parallelType = 'local';
    disp('run in complex1')
  else %OU
    parallelType = 'local';
    disp('run in OU')
  end
end

if ~isfield(auxData,'basicSeqNameDir') 
  auxData.basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];
end

if ~isfield(auxData,'auxData.basicSeqKey')
  auxData.basicSeqKey = [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
end
disp(['working with database: ',auxData.basicSeqNameDir,'; Is it ok?'])



keyboard


currNumProcessors = auxData.numProcessors;

n = length(auxData.dr);
part = 1:currNumProcessors:n;
part(end) = n+1;
tmpInd = cell(length(part)-1,1);
for i=1:length(part)-1
  tmpInd{i} = part(i):part(i+1)-1;
end



if strcmp(parallelType,'local')

  % do not touch this
  if currNumProcessors>1 
   
    
    %w = ['matlabpool open local ',num2str(currNumProcessors),';parfor i=1:length(tmpInd),i,[xx,ss]=solveForGroupOrFourth(i,i,auxData,[],[],tmpInd,basicSeqNameDir,basicSeqKey,readLength,uniqueReads,uniqueReads_length);x{i} = xx{i};sumRelevantReads(i) = ss(i);;end;matlabpool close'];
    %eval(w);
    
  else % single
    %keyboard 
    
    
    for i=1:length(tmpInd)
      for j=tmpInd{i}
        ResultsDistance([userDir,'/CS/BAC/',auxData.outFileDirName,'/'],auxData.dr(j).name,auxData);
      end
      
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
  
  for j=1:length(partDist)-1
    allocateToOne = partDist(j):partDist(j+1)-1;
    mFileName = [runFileName,'_partDist_',num2str(j)];
    saveFileName = [basicAllocationFileName,'_partDist_',num2str(j),'_res'];
    runOneNodeBasedSolveForGroupOrFourth(allocateToOne,basicAllocationFileName,mFileName,saveFileName,userDir,[],auxData.queueName);
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
  
  
  %  Br 
  
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
	
  
  unix(['rm -f ',auxData.currDir,'/core*; ','rm -f ',runFileName,'_*; ','rm -f ',basicSaveName,'_*; ']);
   
  save([basicSaveName,'_tmpRes'],'x','sumRelevantReads')	
	    
end % strcmp(parallelType,'cluster')

