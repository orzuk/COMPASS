function distributeBAC_generalOr(auxDataFileName)

% the same as distributeBAC_general but also does infinite reads. allows differnet normalizations for the number of reads

% reads are red
% can load specific reads - just like needed in similar BAC! just make sure the reads are called red

% create stand alone run files and then run them
%rmpath(genpath('/home/csfaculty/shental/Weizmann/ph2users/fenoam/'));rmpath('/home/csfaculty/shental/matlab');cd ~/CS/mFilesBAC; addpath('~/CS/mFilesBAC'); mcc -m distributeBAC_generalOr  -d /homes/csfaculty/shental/CS/BAC/cvx  -a currReads.m -a findParametersReads.m -a fun_prep_data_for_simulations_loadSpecificBAC2Or.m -a get_duplicates.m -a get_duplicates2.m -a getuserdir.m -a iterateParallelDistributedSequenceFilesOr.m  -a myHistForBAC.m -a myUniqueUINT8.m -a randi.m -a newReadsBasedSeed2Or.m -a currReads.m -a prepareReadsForRandomSequences.m -a solveForGroupOr.m   -a prepareGroupOf1000DistributedSequenceFilesOr  -a runOneGroupOf1000ForCompilation.m -a set_solution_bounds.m
 

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

if ~isfield(auxData,'randomizeSeedForRepeatRandomGroups') % randomizeSeedForRepeatRandomGroups to make each copy different
  auxData.randomizeSeedForRepeatRandomGroups = 0;
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

if ~isfield(auxData,'basicSeqNameDir') 
  auxData.basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];
end
if ~isfield(auxData,'basicSeqKey')
  auxData.basicSeqKey = [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
end
disp(['working with database: ',auxData.basicSeqNameDir,'; Is it ok?'])




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
%keyboard
% create the reads
if auxData.inifiniteNumberOfReadsFlag==0
  disp('finite number of reads')
  if auxData.createReadsFlag 
    [junk,red] = newReadsBasedSeed2Or(auxData.reads_data_fileName,userDir,auxData.basicSeqNameDir,auxData.basicSeqKey,auxData.numBACtoConsider);
  else % load the reads themselves
    load(auxData.reads_data_fileName)
  end
  
  %keyboard
  
  %keyboard
  [uniqueReads,uniqueReads_inds] = extract_sub_kmers(red, auxData.readLength*ones(size(red,1),1),auxData.readLength, 1,0);
  
  [junk_vals junk_inds uniqueReads_length] = get_duplicates(uniqueReads_inds(:,1));
  clear junk_vals junk_inds uniqueReads_inds
     
  
   
else
  disp('using an infinite number of reads');

  load(auxData.reads_data_fileName,'correctWeight','ind_bac_in_mix')
  
  % for infinite
  [uniqueReads,uniqueReads_length,auxData.fracRelevantReadsForInfinity]=createReadsForInfiniteNumberOr(ind_bac_in_mix,correctWeight,auxData.readLength,auxData.basicSeqNameDir,auxData.basicSeqKey);
 
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
      [cX{l},cSumRelevantReads] = iterateParallelDistributedSequenceFilesOr(uniqueReads,uniqueReads_length,curr_kp,auxData.basicSeqNameDir,auxData.basicSeqKey,parallelType,[basicSaveName,'_k_',num2str(k)],userDir,[runFileName,'_k_',num2str(k)],auxData);
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
    disp('running one random split. iterateParallelDistributedSequenceFilesOr.m ')
    [currX,currSumRelevantReads] = iterateParallelDistributedSequenceFilesOr(uniqueReads,uniqueReads_length,curr_kp,auxData.basicSeqNameDir,auxData.basicSeqKey,parallelType,[basicSaveName,'_k_',num2str(k)],userDir,[runFileName,'_k_',num2str(k)],auxData);
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
  [currX,currSumRelevantReads] = iterateParallelDistributedSequenceFilesOr(uniqueReads,uniqueReads_length,curr_kp,auxData.basicSeqNameDir,auxData.basicSeqKey,parallelType,[basicSaveName,'_final_',num2str(i)],userDir,[runFileName,'_final_',num2str(i)],auxData);
  found{i} = zeros(1,auxData.numBACtoConsider);
  found{i}(curr_kp) = currX;
end
save(auxData.saveName,'found','store_kp','store_X','store_sumRelevantReads')


if exist('found')
  unix(['rm -f ',auxData.currDir,'/',auxData.tmpFileName,'_*']);
  %unix(['rm -f ',auxData.currDir,'/run_',auxData.tmpFileName,'_*']);
end


if 1==2
    
  % OU
    
  clear
  load([getuserdir,'/CS/BAC/cellsimbac2'])
  Nbac_in_mixture = 500;
  Nreads = 2e6;
  readlen = 50;
  npower = 0.5;
  bac_dist_flag = 0;
  outfilename = '~/CS/BAC/dataForSim/sim44';
  N = 2001;
  createReadsFlag = 0;
  [strinfilename] = fun_prep_data_for_simulations_loadSpecificBAC2([],[],indsimbac_vec,Nbac_in_mixture,Nreads,readlen,npower,bac_dist_flag,[],outfilename,N,createReadsFlag);

  auxData = struct;
  auxData.inifiniteNumberOfReadsFlag = 0;
  auxData.reads_data_fileName = strinfilename;
  auxData.tmpFileName = 'testNew';
  auxData.currDir = ['/homes/csfaculty/shental/CS/BAC/tmpRuns/',auxData.tmpFileName];
  auxData.basicSeqNameDir = ['~/CS/BAC/datNoNonACGT/packed64/'];
  auxData.basicSeqKey= ['~/CS/BAC/datNoNonACGT/keyNoNonACGT'];
  
  %unix(['mkdir ',auxData.currDir]);
  %unix(['cd ',auxData.currDir]); 
  
  
  
  auxData.readLength = 50;
  auxData.saveName = ['~/CS/BAC/',auxData.tmpFileName];
  auxData.numProcessors = 1;
  auxData.groupSize = 1000;
  auxData.smallestSetOfCollected = 1000;
  auxData.numBACtoConsider = N;
  auxData.thresholdForCollectingBAC = 1e-3;
  auxData.firstFlag = 1;
  auxData.createReadsFlag = 1;
  auxDataFileName = ['~/CS/tmp/structFor_',auxData.tmpFileName];
  save(auxDataFileName,'auxData')
  distributeBAC_generalOr(auxDataFileName)

end


if 1==2 % large
    
  % OU
  profile -memory on;  
  clear
  load([getuserdir,'/CS/BAC/cellsimbac2'])
  Nbac_in_mixture = 500;
  Nreads = 1e5;
  readlen = 50;
  npower = 0.5;
  bac_dist_flag = 0;
  outfilename = '~/CS/BAC/dataForSim/sim44';
  N = 410849;
  createReadsFlag = 0;
  [strinfilename] = fun_prep_data_for_simulations_loadSpecificBAC2([],[],indsimbac_vec,Nbac_in_mixture,Nreads,readlen,npower,bac_dist_flag,[],outfilename,N,createReadsFlag);

  auxData = struct;
  auxData.repeatRandomGroups = 10;
  auxData.moveDependentToNextStageFlag = 1;
  auxData.randomizeSeedAfterCreateReadsFlag = 1;
  
  auxData.inifiniteNumberOfReadsFlag = 0;
  auxData.reads_data_fileName = strinfilename;
  auxData.tmpFileName = 'testNew';
  auxData.currDir = ['/homes/csfaculty/shental/CS/BAC/tmpRuns/',auxData.tmpFileName];
  auxData.basicSeqNameDir = ['~/CS/BAC/datNoNonACGT/packed64/'];
  auxData.basicSeqKey= ['~/CS/BAC/datNoNonACGT/keyNoNonACGT'];
  
  %unix(['mkdir ',auxData.currDir]);
  %unix(['cd ',auxData.currDir]); 
  
  
  
  auxData.readLength = 50;
  auxData.saveName = ['~/CS/BAC/',auxData.tmpFileName];
  auxData.numProcessors = 1;
  auxData.groupSize = 1000;
  auxData.smallestSetOfCollected = 1000;
  auxData.numBACtoConsider = N;
  auxData.thresholdForCollectingBAC = 1e-3;
  auxData.firstFlag = 1;
  auxData.createReadsFlag = 1;
  auxDataFileName = ['~/CS/tmp/structFor_',auxData.tmpFileName];
  save(auxDataFileName,'auxData')
  distributeBAC_generalOr(auxDataFileName);profile report;profile off;
end
