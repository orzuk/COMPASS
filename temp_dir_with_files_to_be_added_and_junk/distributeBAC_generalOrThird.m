function distributeBAC_generalOrThird(auxDataFileName)

% the same as distributeBAC_general but also does infinite reads. allows differnet normalizations for the number of reads

% reads are red
% can load specific reads - just like needed in similar BAC! just make sure the reads are called red

% create stand alone run files and then run them
%rmpath(genpath('/home/csfaculty/shental/Weizmann/ph2users/fenoam/'));rmpath('/home/csfaculty/shental/matlab');cd ~/CS/mFilesBAC;addpath('~/CS/mFilesBAC'); mcc -v -m distributeBAC_generalOrThird  -d /homes/csfaculty/shental/CS/BAC/cvx -a grep.m  -a currReads.m -a findParametersReads.m -a fun_prep_data_for_simulations_loadSpecificBAC2Or.m -a get_duplicates.m -a get_duplicates2.m -a getuserdir.m -a iterateParallelDistributedSequenceFilesOr2.m  -a myHistForBAC.m -a myUniqueUINT8.m -a randi.m -a newReadsBasedSeed2Or.m -a currReads.m -a prepareReadsForRandomSequences.m -a solveForGroupOr.m   -a prepareGroupOf1000DistributedSequenceFilesOr  -a runOneGroupOf1000ForCompilation.m -a set_solution_bounds.m -a doRepeats2.m -a calcPart.m 

if isdeployed  %load a file which has a variable by the name auxData 
  load(auxDataFileName)
else
  load(auxDataFileName)
end


if ~isfield(auxData,'correctMeasurementForReadErrorsFlag')
  auxData.correctMeasurementForReadErrorsFlag = 0;
end

if ~isfield(auxData,'batchSize') % number of runs per run in cluster - timewise!
  auxData.batchSize = 400;
end

if ~isfield(auxData,'upperLimitForRepeat')
  auxData.upperLimitForRepeat = 150000;
end

if ~isfield(auxData,'repeatWhenLowerThanThisValue')
  auxData.repeatWhenLowerThanThisValue = 20000;
end

if ~isfield(auxData,'saveFullDataFlag')
  auxData.saveFullDataFlag = 0;
end

if ~isfield(auxData,'realReadsFromFileFlag')
  auxData.realReadsFromFileFlag = 0;
end

if ~isfield(auxData,'addNoiseFlag')
  auxData.addNoiseFlag = 0;
end

if auxData.addNoiseFlag
  if ~isfield(auxData,'ErrorStruct')
    disp('ErrorStruct??!! distributeBAC_generalOrThird.m')
    return
  end
else
  auxData.ErrorStruct = [];
end

if ~isfield(auxData,'forcePartFromOutsideFlag')
  auxData.forcePartFromOutsideFlag = 0;
end

if ~isfield(auxData,'preprocessByExistanceFlag')
  auxData.preprocessByExistanceFlag = 0;
end
if ~isfield(auxData,'repeatRandomGroups') % repeat random selections of groups 
  auxData.repeatRandomGroups = 1;
end

if ~isfield(auxData,'inifiniteNumberOfReadsFlag')
  auxData.inifiniteNumberOfReadsFlag = 0;
end

if ~isfield(auxData,'keepOriginalOrderFlag') % keepOriginalOrderFlag for  checking in the original order
  auxData.keepOriginalOrderFlag = 0;
end


if ~isfield(auxData,'randomizeSeedAfterCreateReadsFlag') % randomizeSeedAfterCreateReadsFlag for considering a certain bacteria from a group
  auxData.randomizeSeedAfterCreateReadsFlag = 0;
end



% add to auxData stuff
if ~isfield(auxData,'groupSize') % partition into groups of groupSize bacteria
  auxData.groupSize = 1000;
end

if ~isfield(auxData,'smallestSetOfCollected') % after collecting the bacteria from each group - continue iterating id their number is smaller than smallestSetOfCollected
  auxData.smallestSetOfCollected = 1000;
end

if ~isfield(auxData,'thresholdForCollectingBAC') % thresholdForCollectingBAC for considering a certain bacteria from a group
  auxData.thresholdForCollectingBAC = 1e-3;
end

if ~isfield(auxData,'numBACtoConsider')
  auxData.numBACtoConsider = 455055;
end

if ~isfield(auxData,'createReadsFlag')
  auxData.createReadsFlag = 0;
end

%%%%%% %end add to auxData stuff


if ~isfield(auxData,'brFlag')
  auxData.brFlag = 0;
end

if ~isfield(auxData,'queueName')
  auxData.queueName = [];
end

if ~isfield(auxData,'firstFlag')
  auxData.firstFlag = 0;
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
if ~isfield(auxData,'basicSeqKey')
  auxData.basicSeqKey = [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
end
disp(['working with database: ',auxData.basicSeqNameDir,'; Is it ok?'])




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end load variables

%keyboard
fn = fieldnames(auxData);
fprintf('%%%%%%%%%%%\n');
fprintf('data:\n');
for i=1:length(fn)
  if ~strcmp(fn{i},'ErrorStruct')
    if isnumeric(auxData.(fn{i})) 
      fprintf('%s: %d\n',fn{i},auxData.(fn{i}));
    else
      fprintf('%s: %s\n',fn{i},auxData.(fn{i}));
    end
    
  else % ErrorStruct
    fprintf('\n\n');
    if auxData.addNoiseFlag
      fprintf('ErrorStruct was considered\n\n\n')
    end
  end % if ErrorStruct
  
end    
fprintf('end data\n')
fprintf('%%%%%%%%%%%%\n')

%%%%%%%%%%%%%%%%%%555

if ~isdeployed && auxData.firstFlag==1 
  addpath([userDir,'/CS/mFilesBAC/']);
  addpath(genpath([userDir,'/CS/BAC/cvx/']))
end


%keyboard

%unix(['mkdir ',auxData.currDir]);
runFileName = [auxData.currDir,'/run_',auxData.tmpFileName];
basicSaveName = [auxData.currDir,'/',auxData.tmpFileName,'_dat'];
%keyboard

% create the reads or read them from file
if auxData.realReadsFromFileFlag==0
  disp('using simulated reads')
  if auxData.inifiniteNumberOfReadsFlag==0
    disp('finite number of reads')
    %keyboard
    if auxData.createReadsFlag 
      [junk,red] = newReadsBasedSeed2Or(auxData.reads_data_fileName,userDir,auxData.basicSeqNameDir,auxData.basicSeqKey,auxData.numBACtoConsider,auxData.addNoiseFlag,auxData.ErrorStruct);
    else % load the reads themselves
      load(auxData.reads_data_fileName);
    end
    
    %keyboard
    
    %keyboard
    [uniqueReads,uniqueReads_inds] = extract_sub_kmers(red, auxData.readLength*ones(size(red,1),1),auxData.readLength, 1,0);
    
    [junk_vals junk_inds uniqueReads_length] = get_duplicates(uniqueReads_inds(:,1));
    clear junk_vals junk_inds uniqueReads_inds
        
    disp('delete this distributeBAC_generalOrThird.m')
    auxData.readsForCorrectReads = zeros(sum(uniqueReads_length),2,'uint64');
    kk = 0;
    for j=1:length(uniqueReads_length)
      auxData.readsForCorrectReads(kk+1:kk+uniqueReads_length(j),:) = repmat(uniqueReads(j,:),uniqueReads_length(j),1);
      kk = kk+uniqueReads_length(j);
    end
    %%%%%%%%%%%%%%%%%
    
  else
    disp('using an infinite number of reads');
    
    load(auxData.reads_data_fileName,'correctWeight','ind_bac_in_mix')
    
    % for infinite 
    
    [uniqueReads,uniqueReads_length,auxData.fracRelevantReadsForInfinity]=createReadsForInfiniteNumberOr(ind_bac_in_mix,correctWeight,auxData.readLength,auxData.basicSeqNameDir,auxData.basicSeqKey);
    
    clear ind_bac_in_mix correctWeight
  end

else % real reads
  disp('using reads from file')
  load(auxData.reads_data_fileName) % uniqueReads uniqueReads_length
  
end % load reads



if auxData.randomizeSeedAfterCreateReadsFlag
  disp('randomized state, after the reads were created')
  rand('seed',sum(100*clock));
else
  disp('use the original seed')
end

%keyboard
clear store* 
k = 1;

curr_kp = 1:auxData.numBACtoConsider;
while length(curr_kp)>auxData.smallestSetOfCollected 
  %keyboard
  if auxData.repeatRandomGroups>1
    
    if length(curr_kp)>auxData.upperLimitForRepeat
      % regular
      disp('running one random split. iterateParallelDistributedSequenceFilesOr2.m ')
      [currX,currSumRelevantReads] = iterateParallelDistributedSequenceFilesOr2(uniqueReads,uniqueReads_length,curr_kp,auxData.basicSeqNameDir,auxData.basicSeqKey,parallelType,[basicSaveName,'_k_',num2str(k)],userDir,[runFileName,'_k_',num2str(k)],auxData);
    else
      % less than 150000. when lower than 50000 does in one batch
      [currX,currSumRelevantReads] = doRepeats2(uniqueReads,uniqueReads_length,curr_kp,parallelType,basicSaveName,k,userDir,runFileName,auxData);
      %what about the value of k?
    end
    
    
  else % regular
    disp('running one random split. iterateParallelDistributedSequenceFilesOr2.m ')
    [currX,currSumRelevantReads] = iterateParallelDistributedSequenceFilesOr2(uniqueReads,uniqueReads_length,curr_kp,auxData.basicSeqNameDir,auxData.basicSeqKey,parallelType,[basicSaveName,'_k_',num2str(k)],userDir,[runFileName,'_k_',num2str(k)],auxData);
  
  end
  
  
  
  
  next_kp = curr_kp(find(currX>auxData.thresholdForCollectingBAC));
  store_kp{k} = curr_kp;
  store_X{k} = currX;
  store_sumRelevantReads{k} = currSumRelevantReads;
  k = k +1;
  
  curr_kp = next_kp;
  
  if length(store_kp{k-1})==length(curr_kp)
    disp('did not reduce size, so breaks. distributeBAC_general3.m')
    break
  end
end
if length(store_kp{k-1})~=length(curr_kp) % size was reduced
  store_kp{k} = curr_kp;
  store_X{k} = [];
end

disp('does it take the last one ok? distributeBAC_general3.m')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve again
lengths_kp = zeros(1,length(store_kp));
for i=1:length(lengths_kp)
  lengths_kp(i) = length(store_kp{i});
end
%1
%keyboard
do = find(lengths_kp<auxData.smallestSetOfCollected*2);

% deal with the case that n-1 is larger than auxData.smallestSetOfCollected*2 and the n'th is really small
if isempty(find(lengths_kp(do)>100))
  do = do(1)-1;
  disp(['running with larger group than ',num2str(auxData.smallestSetOfCollected*2),' of size: ',num2str(lengths_kp(do)),' distributeBAC_general3.m'])
end


% solve - regular but the negative solution are the ones which are unique_inds
auxData.solveOr = 1;
auxData.do_lp_flag = 0;
auxData.solveFinalGroupFlag = 1;
auxData.numProcessors = 1;
auxData.moveDependentToNextStageFlag = 0;
%auxData.keepOriginalOrderFlag==1;
for i=do
  curr_kp = store_kp{i};
  auxData.groupSize = length(curr_kp)-1;
  [currX,currSumRelevantReads] = iterateParallelDistributedSequenceFilesOr2(uniqueReads,uniqueReads_length,curr_kp,auxData.basicSeqNameDir,auxData.basicSeqKey,parallelType,[basicSaveName,'_final_',num2str(i)],userDir,[runFileName,'_final_',num2str(i)],auxData);
  found{i} = zeros(1,auxData.numBACtoConsider);
  found{i}(curr_kp) = currX;
end

if auxData.saveFullDataFlag
  load(auxData.reads_data_fileName,'correctWeight','ind_bac_in_mix')  
  disp('standard way of saving')
  save(auxData.saveName,'found','store_kp','store_X','store_sumRelevantReads','correctWeight','ind_bac_in_mix')
else
  disp('save compact way with the correctWeight')
%keyboard
  load(auxData.reads_data_fileName,'correctWeight','ind_bac_in_mix')  
  nonEmpty = find(cellfun(@isempty,found)==0);
  resCell = cell(1+length(nonEmpty),1);
  resCell{1} = [ind_bac_in_mix,correctWeight(ind_bac_in_mix)];
  for zz=1:length(nonEmpty)
    nonZer = find(abs(found{nonEmpty(zz)})>10^-5);
    resCell{zz+1} = [nonZer',found{nonEmpty(zz)}(nonZer)'];
  end
  save_store_kp = store_kp(2:end);
  save_store_X = store_X(2:end);
  save(auxData.saveName,'resCell','save_store_kp','save_store_X');
end

%keyboard
if exist('found')
  %unix(['rm -f ',auxData.currDir,'/',auxData.tmpFileName,'_*']);
  unix(['rm -fr ',auxData.currDir,'/']);
  %unix(['rm -f ',auxData.currDir,'/run_',auxData.tmpFileName,'_*']);
end


if 1==2 % large
  clear
  userDir = getuserdir;
  outFileDirName = 'ffAddPrimer';
  basicName = outFileDirName;
  outfilename = [userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',basicName];
  unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName,...
       ';mkdir ',userDir,'/CS/BAC/',outFileDirName])
  
  load([userDir,'/CS/BAC/cellsimbac2'])
  Nbac_in_mixture = 100;
  Nreads = Inf;
  readlen = 400;
  npower = 0;
  bac_dist_flag = 0;
  
  N = 410849;
  createReadsFlag = 0;
  [strinfilename] = fun_prep_data_for_simulations_loadSpecificBAC2([],[],indsimbac_vec,Nbac_in_mixture,Nreads,readlen,npower,bac_dist_flag,[],outfilename,N,createReadsFlag);
  
  
  
  
  auxData = struct;

  auxData.repeatRandomGroups = 10;
  auxData.moveDependentToNextStageFlag = 1;
  
  auxData.tmpFileName = outFileDirName;
  auxData.randomizeSeedAfterCreateReadsFlag = 1;
  auxData.saveName = [userDir,'/CS/BAC/',outFileDirName,'/',auxData.tmpFileName];
  auxData.inifiniteNumberOfReadsFlag = 1;
  auxData.currDir = [userDir,'/CS/BAC/tmpRuns/',outFileDirName,'/',auxData.tmpFileName];
  auxData.basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];
  auxData.basicSeqKey= [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
  auxData.createReadsFlag = 1;
  auxData.brFlag =  0;
  auxData.queueName = '';  
  if auxData.brFlag
    auxData.reads_data_fileName = strinfilename;
  else
    tmpStr = strrep(strinfilename,userDir,userDir);
    auxData.reads_data_fileName = tmpStr;
  end
  
  auxData.readLength = readlen;
  auxData.numProcessors = 6;
  auxData.firstFlag = 0;
  auxData.groupSize = 1000;
  auxData.smallestSetOfCollected = 1000;
  auxData.numBACtoConsider = N;
  auxData.thresholdForCollectingBAC = 1e-3;
  auxDataFileName = [userDir,'/CS/tmp/',outFileDirName,'/','structFor_',auxData.tmpFileName];
  auxData.addNoiseFlag = 0;
  auxData.correctMeasurementForReadErrorsFlag = 0;
  auxData.saveFullDataFlag = 1;
  auxData.upperLimitForRepeat = 150000;
  auxData.repeatWhenLowerThanThisValue = 20000;
  auxData.batchSize = 400;
  
  save(auxDataFileName,'auxData')
  
  distributeBAC_generalOrThird(auxDataFileName);
  
end

if 1==2 %small
  clear
  userDir = getuserdir;
  outFileDirName = 'ffSmall';
  basicName = outFileDirName;
  outfilename = [userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',basicName];
  unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,...
        ';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,...
        ';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName])
  
  load([userDir,'/CS/BAC/cellsimbac2'])
  Nbac_in_mixture = 10;
  Nreads = Inf;
  readlen = 50;
  npower = 0;
  bac_dist_flag = 0;
  
  N = 10000;
  createReadsFlag = 0;
  [strinfilename] = fun_prep_data_for_simulations_loadSpecificBAC2([],[],indsimbac_vec,Nbac_in_mixture,Nreads,readlen,npower,bac_dist_flag,[],outfilename,N,createReadsFlag);
  
  
  
  
  auxData = struct;

  auxData.repeatRandomGroups = 10;
  auxData.moveDependentToNextStageFlag = 1;
  
  auxData.tmpFileName = outFileDirName;
  auxData.randomizeSeedAfterCreateReadsFlag = 1;
  auxData.saveName = [userDir,'/CS/BAC/',outFileDirName,'/',auxData.tmpFileName];
  auxData.inifiniteNumberOfReadsFlag = 1;
  auxData.currDir = [userDir,'/CS/BAC/tmpRuns/',outFileDirName,'/',auxData.tmpFileName];
  auxData.basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];
  auxData.basicSeqKey= [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
  auxData.createReadsFlag = 1;
  auxData.brFlag =  0;
  auxData.queueName = '';  
  if auxData.brFlag
    auxData.reads_data_fileName = strinfilename;
  else
    tmpStr = strrep(strinfilename,userDir,userDir);
    auxData.reads_data_fileName = tmpStr;
  end

  auxData.readLength = 50;
  auxData.numProcessors = 7;
  auxData.firstFlag = 0;
  auxData.groupSize = 1000;
  auxData.smallestSetOfCollected = 100;
  auxData.numBACtoConsider = N;
  auxData.thresholdForCollectingBAC = 1e-3;
  auxDataFileName = [userDir,'/CS/tmp/',outFileDirName,'/','structFor_',auxData.tmpFileName];
  
  auxData.addNoiseFlag = 0;
  auxData.correctMeasurementForReadErrorsFlag = 0;
  auxData.saveFullDataFlag = 1;
  auxData.upperLimitForRepeat = 6000;
  auxData.repeatWhenLowerThanThisValue = 2000;
  auxData.batchSize = 20;
  save(auxDataFileName,'auxData')
  
  distributeBAC_generalOrThird(auxDataFileName);
  
end

