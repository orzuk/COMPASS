function distributeBAC_general3(auxDataFileName)

% the same as distributeBAC_general but also does infinite reads. allows differnet normalizations for the number of reads

% reads are red
% can load specific reads - just like needed in similar BAC! just make sure the reads are called red

% create stand alone run files and then run them
%rmpath(genpath('/home/csfaculty/shental/Weizmann/ph2users/fenoam/'));rmpath('/home/csfaculty/shental/matlab');cd ~/CS/mFilesBAC; addpath('~/CS/mFilesBAC'); mcc -m distributeBAC_general3  -d /homes/csfaculty/shental/CS/BAC/cvx  -a currReads.m -a findParametersReads.m -a fun_prep_data_for_simulations_loadSpecificBAC2.m -a get_duplicates.m -a get_duplicates2.m -a getuserdir.m -a iterateParallelDistributedSequenceFiles3.m  -a myHistForBAC.m -a myUniqueUINT8.m -a randi.m -a newReadsBasedSeed2.m -a currReads.m -a prepareReadsForRandomSequences.m -a solveForGroup3.m   -a prepareGroupOf1000DistributedSequenceFiles2  -a runOneGroupOf1000ForCompilation.m -a set_solution_bounds.m

%cd ~/CS; addpath('~/CS/mFilesBAC'); mcc -m distributeBAC_general4  -d /ph2users/fenoam/CS/BAC/cvx  

%keyboard

%cd /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/tmpRuns/testSimSec500MoveDependentFinite_round1001 ; /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/cvx/run_distributeBAC_general3.sh /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/structFor_testSimSec500MoveDependentFinite_round1001 ;

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
while length(curr_kp)>auxData.smallestSetOfCollected 
  %keyboard
  if length(curr_kp)<20000 & auxData.repeatRandomGroups>1
    disp(['running with: ',num2str(length(curr_kp)),'; repeating runs: ',num2str(auxData.repeatRandomGroups),' distributeBAC_general3.m'])
    n_kp = zeros(auxData.repeatRandomGroups,length(curr_kp));
    mx = zeros(1,length(curr_kp));
    mx_currSumRelevantReads = 0;
    cSumRelevantReads = 0;
    for l=1:auxData.repeatRandomGroups
      [cX{l},cSumRelevantReads] = iterateParallelDistributedSequenceFiles3(uniqueReads,uniqueReads_length,curr_kp,basicSeqNameDir,basicSeqKey,parallelType,[basicSaveName,'_k_',num2str(k)],userDir,[runFileName,'_k_',num2str(k)],auxData);
      n_kp(l,find(cX{l}>auxData.thresholdForCollectingBAC)) = 1;
      
      mx = max(mx,cX{l});
      mx_currSumRelevantReads = max(mx_currSumRelevantReads,cSumRelevantReads);
    end
    
    sum_n_kp = sum(n_kp,1);
    keepMajority = find(sum_n_kp>=auxData.repeatRandomGroups/2);

    currX = zeros(1,length(curr_kp));
    currX(keepMajority) = mx(keepMajority);
    currSumRelevantReads = cSumRelevantReads;
    
  else
    disp('running one random split. iterateParallelDistributedSequenceFiles3.m ')
    [currX,currSumRelevantReads] = iterateParallelDistributedSequenceFiles3(uniqueReads,uniqueReads_length,curr_kp,basicSeqNameDir,basicSeqKey,parallelType,[basicSaveName,'_k_',num2str(k)],userDir,[runFileName,'_k_',num2str(k)],auxData);
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
  [currX,currSumRelevantReads] = iterateParallelDistributedSequenceFiles3(uniqueReads,uniqueReads_length,curr_kp,basicSeqNameDir,basicSeqKey,parallelType,[basicSaveName,'_final_',num2str(i)],userDir,[runFileName,'_final_',num2str(i)],auxData);
  found{i} = zeros(1,auxData.numBACtoConsider);
  found{i}(curr_kp) = currX;
end
save(auxData.saveName,'found','store_kp','store_X','store_sumRelevantReads')


if exist('found')
  unix(['rm -f ',auxData.currDir,'/',auxData.tmpFileName,'_*']);
  %unix(['rm -f ',auxData.currDir,'/run_',auxData.tmpFileName,'_*']);
end
