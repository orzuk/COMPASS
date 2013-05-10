function [X,sumRelevantReads]=iterateParallelDistributedSequenceFiles2(uniqueReads,uniqueReads_length,curr_kp,basicSeqNameDir,basicSeqKey,parallelType,basicSaveName,userDir,runFileName,auxData)
%keyboard
% deleted the following input as they appear in auxData: readLength,numProcessors,,queueName

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

disp('should change the cluste according to the changes in local')
if strcmp(parallelType,'local')
  %keyboard
  %%%%%%%%
  % do not touch this
  if currNumProcessors>1 
    %keyboard
    w = ['matlabpool open local ',num2str(currNumProcessors),';parfor i=1:length(tmpInd),i,[normalizedBac values] = prepareGroupOf1000DistributedSequenceFiles2(readLength,tmpInd{i},basicSeqNameDir,basicSeqKey);[fracRelevantReads,sumRelevantReads(i)] = currReads(uniqueReads,uniqueReads_length,values);if sumRelevantReads(i)>0,[x{i},firstFlag]=runOneGroupOf1000post(1,normalizedBac,fracRelevantReads,0);else,x{i} = zeros(size(normalizedBac,2),1);end;end;matlabpool close'];
    eval(w);
    
  
  else % single
    %keyboard 
    for i=1:length(tmpInd)
       i
       % create a part of the matrix
       [normalizedBac values] = prepareGroupOf1000DistributedSequenceFiles2(readLength,tmpInd{i},basicSeqNameDir,basicSeqKey);  
       
       [fracRelevantReads,sumRelevantReads(i)] = currReads(uniqueReads,uniqueReads_length,values);
       clear values
       if sumRelevantReads(i)>0
         % find for 1000, output the sum of relevant reads 
         [x{i},firstFlag]=runOneGroupOf1000post(1,normalizedBac,fracRelevantReads,0);
       else
         x{i} = zeros(size(normalizedBac,2),1);
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
  
  
  % create the reads locations and sumRelevantReads
  basicAllocationFileName = [basicSaveName,'_basicAllocationData'];
  save(basicAllocationFileName,'readLength','basicSeqNameDir','basicSeqKey','uniqueReads','uniqueReads_length','tmpInd')
  
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
    if 1==2  % old way of submitting jobs
      for j=1:length(partDist)-1
        if j>=0
          pause(30)
          disp('pause for 30 seconds')
        end
        
        allocateToOne = partDist(j):partDist(j+1)-1;
        mFileName = [basicAllocationFileName,'_partDist_',num2str(j)];
        saveFileName = [basicAllocationFileName,'_partDist_',num2str(j),'_res'];
        saveFileName
        runOneNodeCompiled2(allocateToOne,basicAllocationFileName,mFileName,saveFileName,userDir,[],auxData.queueName);
      end
    end
    
	%keyboard
    % new way of submitting jobs - create the runFiles
    for j=1:length(partDist)-1
      allocateToOne = partDist(j):partDist(j+1)-1;
      mFileName = [runFileName,'_partDist_',num2str(j)];
      saveFileName = [basicAllocationFileName,'_partDist_',num2str(j),'_res'];
      runOneNodeCompiled3(allocateToOne,basicAllocationFileName,mFileName,saveFileName,userDir,[],auxData.queueName);
    end
    %keyboard
    % chmod
    unix(['chmod 0700 ',runFileName,'_partDist_*']);
    % run them
    
    tmpRunFileName = runFileName;
    delSlash = find(tmpRunFileName=='/');
    tmpRunFileName = tmpRunFileName(delSlash(end)+1:end);
    %keyboard
    
    fid = fopen([runFileName,'_shell.sh'],'w');
    fprintf(fid,'#!/bin/bash\n');
    fprintf(fid,'TMPDIR=/tmp/${RANDOM}_${RANDOM}\n');
    fprintf(fid,'mkdir $TMPDIR\n');
    fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
    fprintf(fid,'bsub   -q  %s -J "%s[1-%d]" -o "%s_partDist_Z%s.o" ~/runInput_donotdelete %s_partDist_%s;\n',auxData.queueName,tmpRunFileName,length(partDist)-1,runFileName,'%I',runFileName,'\$LSB_JOBINDEX');
    fprintf(fid,'rm -rf $TMPDIR \n');
    fclose(fid)
    %keyboard
    unix(['chmod 0700 ',runFileName,'_shell.sh;',runFileName,'_shell.sh'])
    
    %unix(['bsub   -q  ',auxData.queueName,...
    %      ' -J "',tmpRunFileName,'[1-',num2str(length(partDist)-1),']"',...
    %      ' -o "',runFileName,'_partDist_Z','%I.o"',...
    %      ' ~/runInput_donotdelete ',...
    %     runFileName,'_partDist_\$LSB_JOBINDEX;']);
   
     %     set myRandomNumber1=`awk ''BEGIN {srand();print int(rand()*50000)}''`;',';set TMPDIR=/tmp/${myRandomNumber1}_${myRandomNumber1};mkdir $TMPDIR;set MCR_CACHE_ROOT=$TMPDIR;...
          
    % bsub   -q hour   -J "myJob[1-3]" -o "/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/run_partDist_%I.o" ./tmpRun.sh  /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/run_partDist_\$LSB_JOBINDEX 

    % check whether it got thrown away
    numSuccessfullyCompleted = zeros(length(partDist)-1,1);
    while length(find(numSuccessfullyCompleted))<length(partDist)-1
      % find succesful - although they may have not run
      unix(['grep -l "Successfully completed" ',runFileName,'_partDist_Z*.o | cut -f2 -d"Z" | cut -f1 -d"." >',runFileName,'_logOfSuccessful.txt;sleep   2'])
      listOfsuccessful = load([runFileName,'_logOfSuccessful.txt']);
      
      
      
      % find unsuccessful and run them again - either did not run or exited
      
      % exit code
      unix(['grep -l "exit code" ',runFileName,'_partDist_Z*.o | cut -f2 -d"Z" | cut -f1 -d"." >',runFileName,'_logOfUnsuccessful1.txt;sleep 2'])
      listOfProblematicBasedExitCode1 = load([runFileName,'_logOfUnsuccessful1.txt']);
      listOfProblematicBasedExitCode1 = listOfProblematicBasedExitCode1';
      
       unix(['grep -l "fatalexit" ',runFileName,'_partDist_Z*.o | cut -f2 -d"Z" | cut -f1 -d"." >',runFileName,'_logOfUnsuccessful2.txt;sleep 2'])
      listOfProblematicBasedExitCode2 = load([runFileName,'_logOfUnsuccessful2.txt']);
      listOfProblematicBasedExitCode2 = listOfProblematicBasedExitCode2';
      
      listOfProblematicBasedExitCode = [listOfProblematicBasedExitCode1;listOfProblematicBasedExitCode2];
      
      % did not run
      listOfProblematicDidNotRun = [];
      for j=listOfsuccessful'
        dr = dir([basicAllocationFileName,'_partDist_',num2str(j),'_res.mat']);
        if length(dr)==1
          numSuccessfullyCompleted(j) = 1;
        else
          listOfProblematicDidNotRun = [listOfProblematicDidNotRun,j];
        end
        

      end
      
      listOfProblematic = [listOfProblematicBasedExitCode;listOfProblematicDidNotRun'];
      
      if ~isempty(listOfProblematic)
        pause(2)
        % delete the .o files and create names
        delThese = [];
        runThese = [];
        for j=1:length(listOfProblematic)
          delThese = [delThese,' ',runFileName,'_partDist_Z',num2str(listOfProblematic(j)),'.o'];
          runThese = [runThese,' ',num2str(listOfProblematic(j))];
        end
        %runThese(1) = [];
        unix(['rm -f ',delThese]);
      
        % submit them again
        tmpName = num2str(runThese);
        tmpName(find(tmpName==' ')) = '_';
        fid = fopen([runFileName,tmpName,'_shell.sh'],'w');
        fprintf(fid,'#!/bin/bash\n');
        fprintf(fid,'TMPDIR=/tmp/${RANDOM}_${RANDOM}\n');
        fprintf(fid,'mkdir $TMPDIR\n');
        fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
        fprintf(fid,'for i in %s; do bsub   -q %s -o %s_partDist_Z%s.o ~/runInput_donotdelete %s_partDist_%s;sleep 10;done;',runThese,auxData.queueName,runFileName,'$i',runFileName,'$i');
        %fprintf(fid,'bsub   -q  %s -J "%s[%s]" -o "%s_partDist_Z%s.o" ~/runInput_donotdelete %s_partDist_%s;\n',auxData.queueName,tmpRunFileName,num2str(runThese),runFileName,'%I',runFileName,'\$LSB_JOBINDEX');
        fprintf(fid,'rm -rf $TMPDIR \n');
        fclose(fid)
        pause(10)
        unix(['chmod 0700 ',runFileName,tmpName,'_shell.sh']);
        disp('does it run?')
        pause(5)
        unix([runFileName,tmpName,'_shell.sh'])
        pause(60)
        
        
        
        
      else % wait for results/problems
        pause(2)
      end
        
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

X = zeros(1,n);
for i=1:length(part)-1
  tmpInd = part(i):part(i+1)-1;
  tmpInd = currInd(tmpInd);
  if ~isempty(x{i})
    X(tmpInd) = x{i};
  end
end




  
