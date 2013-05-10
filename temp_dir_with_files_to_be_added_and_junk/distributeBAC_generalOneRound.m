function distributeBAC_generalOneRound(auxDataFileName)

% the same as distributeBAC_general but also does infinite reads. allows differnet normalizations for the number of reads

% reads are red
% can load specific reads - just like needed in similar BAC! just make sure the reads are called red

% create stand alone run files and then run them
%rmpath(genpath('/home/csfaculty/shental/Weizmann/ph2users/fenoam/'));rmpath('/home/csfaculty/shental/matlab');cd ~/CS/mFilesBAC; addpath('~/CS/mFilesBAC'); mcc -m distributeBAC_generalOneRound  -d /homes/csfaculty/shental/CS/BAC/cvx  -a currReads.m -a findParametersReads.m -a fun_prep_data_for_simulations_loadSpecificBAC2.m -a get_duplicates.m -a get_duplicates2.m -a getuserdir.m -a iterateParallelDistributedSequenceFiles3.m  -a myHistForBAC.m -a myUniqueUINT8.m -a randi.m -a newReadsBasedSeed2.m -a currReads.m -a prepareReadsForRandomSequences.m -a solveForGroup3.m   -a prepareGroupOf1000DistributedSequenceFiles2  -a runOneGroupOf1000ForCompilation.m -a set_solution_bounds.m

%cd ~/CS; addpath('~/CS/mFilesBAC'); mcc -m distributeBAC_general4  -d /ph2users/fenoam/CS/BAC/cvx  

%keyboard

if isdeployed  %load a file which has a variable by the name auxData 
  load(auxDataFileName)
else
  load(auxDataFileName)
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
  if ~isempty(findstr(userDir,'fenoam')) %WIS
    parallelType = 'cluster';
    disp('runs in WIS')
  else %OU
    parallelType = 'local';
    disp('run in OU')
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end load variables

fn = fieldnames(auxData);
fprintf('%%%%%%%%%%%\n');
fprintf('data:\n');
for i=1:length(fn)
  if isnumeric(auxData.(fn{i}))
    fprintf('%s: %d\n',fn{i},auxData.(fn{i}));
  else
    fprintf('%s: %s\n',fn{i},auxData.(fn{i}));
  end
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
basicSeqNameDir = [userDir,'/CS/BAC/dat450000/'];
basicSeqKey =  [userDir,'/CS/BAC/dat450000/key450000'];
%keyboard
% create the reads
if auxData.inifiniteNumberOfReadsFlag==0
  disp('finite number of reads')
  if auxData.createReadsFlag 
    [junk,red] = newReadsBasedSeed2(auxData.reads_data_fileName,userDir,basicSeqNameDir,basicSeqKey,auxData.numBACtoConsider);
  else % load the reads themselves
    load(auxData.reads_data_fileName)
  end
  disp('assumes that red is uint8');
  
  if ~strcmp(class(red),'uint8')  % Check if data is a uint8
    red = uint8(red);         
    disp('converted red to uint8- should have been saved are uint8!!! in distributeBAC_general4.m')
  end
  
  %keyboard
  [uniqueReads,uniqueReads_inds] = myUniqueUINT8(red);
  clear red
  
  uniqueReads_length = zeros(size(uniqueReads,1),1);
  for i=1:size(uniqueReads_length,1)
    uniqueReads_length(i) = length(uniqueReads_inds{i});
  end
  clear uniqueReads_inds
  
else
  disp('using an infinite number of reads');

  load(auxData.reads_data_fileName,'correctWeight','ind_bac_in_mix')
  
  % for infinite
  [uniqueReads,uniqueReads_length,auxData.fracRelevantReadsForInfinity]=createReadsForInfiniteNumber(ind_bac_in_mix,correctWeight,auxData.readLength,basicSeqNameDir,basicSeqKey);
 
  clear ind_bac_in_mix correctWeight
end



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
if auxData.repeatRandomGroups>1
  disp('problem. distributeBAC_generalOneRound')
  keyboard
end

if auxData.preprocessByExistanceFlag==0
  disp('problem. should be preprocessByExistanceFlag=1. distributeBAC_generalOneRound.m')
  keyboard
end

[currX,currSumRelevantReads] = iterateParallelDistributedSequenceFiles3(uniqueReads,uniqueReads_length,curr_kp,basicSeqNameDir,basicSeqKey,parallelType,[basicSaveName,'_k_',num2str(k)],userDir,[runFileName,'_k_',num2str(k)],auxData);
%keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve again
auxData.preprocessByExistanceFlag = 0;

save(auxData.saveName,'currX')

if exist('currX')
  unix(['rm -f ',auxData.currDir,'/',auxData.tmpFileName,'_*']);
  %unix(['rm -f ',auxData.currDir,'/run_',auxData.tmpFileName,'_*']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

