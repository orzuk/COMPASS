function distributeBAC_generalOrSecond(auxDataFileName)

% the same as distributeBAC_general but also does infinite reads. allows differnet normalizations for the number of reads

% reads are red
% can load specific reads - just like needed in similar BAC! just make sure the reads are called red

% create stand alone run files and then run them
%rmpath(genpath('/home/csfaculty/shental/Weizmann/ph2users/fenoam/'));rmpath('/home/csfaculty/shental/matlab');cd ~/CS/mFilesBAC;addpath('~/CS/mFilesBAC'); mcc -v -m distributeBAC_generalOrSecond  -d /homes/csfaculty/shental/CS/BAC/cvx -a grep.m  -a currReads.m -a findParametersReads.m -a fun_prep_data_for_simulations_loadSpecificBAC2Or.m -a get_duplicates.m -a get_duplicates2.m -a getuserdir.m -a iterateParallelDistributedSequenceFilesOr2.m  -a myHistForBAC.m -a myUniqueUINT8.m -a randi.m -a newReadsBasedSeed2Or.m -a currReads.m -a prepareReadsForRandomSequences.m -a solveForGroupOr.m   -a prepareGroupOf1000DistributedSequenceFilesOr  -a runOneGroupOf1000ForCompilation.m -a set_solution_bounds.m
 

if isdeployed  %load a file which has a variable by the name auxData 
  load(auxDataFileName)
else
  load(auxDataFileName)
end

if ~isfield(auxData,'correctMeasurementForReadErrorsFlag')
  auxData.correctMeasurementForReadErrorsFlag = 0;
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
    disp('ErrorStruct??!! distributeBAC_generalOrSecond.m')
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
        
    disp('delete this distributeBAC_generalOrSecond.m')
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
  if length(curr_kp)<auxData.repeatWhenLowerThanThisValue & auxData.repeatRandomGroups>1
    disp(['running with: ',num2str(length(curr_kp)),'; repeating runs: ',num2str(auxData.repeatRandomGroups),' distributeBAC_general3.m'])
      % assumes curr_kp is sorted
    
    partTemplate = 1:auxData.groupSize:length(curr_kp);
    partTemplate(end) = length(curr_kp)+1;
    
    curr_perm = randperm(length(curr_kp));
    curr_kp_forRepeatRandomGroups = curr_kp(curr_perm);
    currPart = partTemplate;
    
    %keyboard
    z = zeros(length(curr_kp),auxData.repeatRandomGroups);
    z((curr_perm),1) = 1:length(curr_kp);
   
    
    for l=2:auxData.repeatRandomGroups
      curr_perm = randperm(length(curr_kp));
      curr_kp_forRepeatRandomGroups = [curr_kp_forRepeatRandomGroups,curr_kp(curr_perm)];
      z((curr_perm),l) = 1+(l-1)*length(curr_kp):length(curr_kp)+(l-1)*length(curr_kp);
      currPart = [currPart,max(currPart)+partTemplate(2:end)-1];
    end
    
    % former version
    %[junk,ind_forRepeatRandoGroups] = get_duplicates(curr_kp_forRepeatRandomGroups);
    %if length(junk)~=length(curr_kp) || ~isempty(find(junk-curr_kp))
    %  disp('problem. repeatRandomGroups in distributeBAC_generalOrSecond.m')
    %  keyboard
    %end
    %z = cell2mat(ind_forRepeatRandoGroups);
    
       
    auxData.keepOriginalOrderFlag = 1;
    auxData.forcePartFromOutsideFlag = 1;
    auxData.partFromOutside = currPart;
    
    [cX,cSumRelevantReads] = iterateParallelDistributedSequenceFilesOr2(uniqueReads,uniqueReads_length,curr_kp_forRepeatRandomGroups,auxData.basicSeqNameDir,auxData.basicSeqKey,parallelType,[basicSaveName,'_k_',num2str(k)],userDir,[runFileName,'_k_',num2str(k)],auxData);
    
    auxData.keepOriginalOrderFlag = 0;
    auxData.forcePartFromOutsideFlag = 0;
    auxData = rmfield(auxData,'partFromOutside');
   
    %keyboard
    value_forRepeatRandomGroups = cX(z);
    if 1==2
      majority_forRepeatRandomGroups = zeros(size(value_forRepeatRandomGroups),'uint16');
      majority_forRepeatRandomGroups(find(value_forRepeatRandomGroups>auxData.thresholdForCollectingBAC)) = 1;
      
      mx_majority = max(value_forRepeatRandomGroups,[],2);
      sum_majority = sum(majority_forRepeatRandomGroups,2);
      keepMajority = find(sum_majority>=auxData.repeatRandomGroups/2);
    else
      disp('average the repeatRandomGroups in distributeBAC_generalOrSecond.m. Does it work well?')
      mx_majority = mean(value_forRepeatRandomGroups,2);
      keepMajority = find(mx_majority>auxData.thresholdForCollectingBAC);
    end
    
    
    currX = zeros(1,length(curr_kp));
    currX(keepMajority) = mx_majority(keepMajority);
    currSumRelevantReads = 0;
    
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
  
  save(auxData.saveName,'resCell');
end

%keyboard
if exist('found')
  %unix(['rm -f ',auxData.currDir,'/',auxData.tmpFileName,'_*']);
  unix(['rm -fr ',auxData.currDir,'/']);
  %unix(['rm -f ',auxData.currDir,'/run_',auxData.tmpFileName,'_*']);
end


if 1==2
    
  % OU
  profile -memory on
  addpath('~/CS/mFilesBAC')
  userDir = getuserdir;
 
  clear
  load([getuserdir,'/CS/BAC/cellsimbac2'])
  Nbac_in_mixture = 1000;
  Nreads = 2e5;
  readlen = 50;
  npower = 0.5;
  bac_dist_flag = 0;
  outfilename = '~/CS/BAC/dataForSim/sim44';
  N = 2*10^4;
  createReadsFlag = 0;
  [strinfilename] = fun_prep_data_for_simulations_loadSpecificBAC2([],[],indsimbac_vec,Nbac_in_mixture,Nreads,readlen,npower,bac_dist_flag,[],outfilename,N,createReadsFlag);

  auxData = struct;
  auxData.inifiniteNumberOfReadsFlag = 1;
  auxData.reads_data_fileName = strinfilename;
  auxData.tmpFileName = 'testNew';
  auxData.currDir = ['/homes/csfaculty/shental/CS/BAC/tmpRuns/',auxData.tmpFileName];
  auxData.basicSeqNameDir = ['~/CS/BAC/datNoNonACGT/packed64/'];
  auxData.basicSeqKey= ['~/CS/BAC/datNoNonACGT/keyNoNonACGT'];
  auxData.repeatRandomGroups = 2;
  
  
  auxData.readLength = 50;
  auxData.saveName = ['~/CS/BAC/',auxData.tmpFileName];
  auxData.numProcessors = 1;
  auxData.groupSize = 1000;
  auxData.smallestSetOfCollected = 1000;
  auxData.numBACtoConsider = N;
  auxData.thresholdForCollectingBAC = 1e-3;
  auxData.firstFlag = 1;
  auxData.createReadsFlag = 1;
  
  % add noise
  
  auxData.ErrorStruct = ErrorStruct;
  
  auxData.correctMeasurementForReadErrorsFlag = 1; % apply Hamming of distance 1
  
  
  
  
  auxDataFileName = ['~/CS/tmp/structFor_',auxData.tmpFileName];
  save(auxDataFileName,'auxData')
  distributeBAC_generalOrSecond(auxDataFileName)

end

if 1==2 % real data
    
  addpath('~/CS/mFilesBAC')
  userDir = getuserdir;
  % OU
  clear
  outfilename = '~/CS/BAC/dataForSim/data43';
  N = 410849;
  createReadsFlag = 0;
  
  auxData = struct;
  auxData.repeatRandomGroups = 10;
  auxData.moveDependentToNextStageFlag = 1;
  auxData.randomizeSeedAfterCreateReadsFlag = 1;
  auxData.realReadsFromFileFlag = 1;
  %auxData.inifiniteNumberOfReadsFlag = 0;
  auxData.reads_data_fileName = '~/SRX059803/output/SRR19843_data';
  auxData.tmpFileName = 'testNewtest43';
  auxData.currDir = ['/homes/csfaculty/shental/CS/BAC/tmpRuns/',auxData.tmpFileName];
  auxData.basicSeqNameDir = ['~/CS/BAC/datNoNonACGT/packed64/'];
  auxData.basicSeqKey= ['~/CS/BAC/datNoNonACGT/keyNoNonACGT'];
  
  %unix(['mkdir ',auxData.currDir]);
  %unix(['cd ',auxData.currDir]); 
  
  
  
  auxData.readLength = 76;
  auxData.saveName = ['~/CS/BAC/',auxData.tmpFileName];
  auxData.numProcessors = 7;
  auxData.groupSize = 1000;
  auxData.smallestSetOfCollected = 1000;
  auxData.numBACtoConsider = N;
  auxData.thresholdForCollectingBAC = 1e-3;
  auxData.firstFlag = 1;
  auxData.createReadsFlag = 1;
  auxDataFileName = ['~/CS/tmp/structFor_',auxData.tmpFileName];
  save(auxDataFileName,'auxData')
  distributeBAC_generalOrSecond(auxDataFileName);profile report;profile off;

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
  distributeBAC_generalOrSecond(auxDataFileName);profile report;profile off;

end

if 1==2
  auxDataFileName = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/tVarNumOfBac/structFor_s8_bac_dist_flag_0_Nbac_in_mixture_300_npower_0_readlen_50_numIter_8_Nreads_100000';
  
  profile -memory on
  distributeBAC_generalOrSecond(auxDataFileName);;profile report;profile off;;profile report;profile off;
  
end

if 1==2
  clear
  load ~/tmp/structFor_ff2_bac_dist_flag_0_Nbac_in_mixture_10_npower_0_readlen_50_numIter_42_Nreads_Inf_Noise_0_Correction_0
  load ~/tmp/ff2_bac_dist_flag_0_Nbac_in_mixture_10_npower_0_readlen_50_numIter_42_Nreads_Inf_Noise_0_Correction_0_Nbacmix_10_Nread_Inf_Readlen_50_Npower_0_bacdistflag_0_NoReads

  userdir = getuserdir;
  auxData.saveName = [userdir,'/CS/BAC/ff2/ff2_bac_dist_flag_0_Nbac_in_mixture_10_npower_0_readlen_50_numIter_42_Nreads_Inf_Noise_0_Correction_0'];
  auxData.saveName = [userdir,'/CS/BAC/ff2/ff2_bac_dist_flag_0_Nbac_in_mixture_10_npower_0_readlen_50_numIter_42_Nreads_Inf_Noise_0_Correction_0'];
  auxData.currDir = [userdir,'/CS/BAC/tmpRuns/ff2/ff2_bac_dist_flag_0_Nbac_in_mixture_10_npower_0_readlen_50_numIter_42_Nreads_Inf_Noise_0_Correction_0'];
  auxData.basicSeqNameDir = [userdir,'/CS/BAC/datNoNonACGT/packed64/'];
  auxData.basicSeqKey = [userdir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
  auxData.reads_data_fileName = [userdir,'/tmp/ff2_bac_dist_flag_0_Nbac_in_mixture_10_npower_0_readlen_50_numIter_42_Nreads_Inf_Noise_0_Correction_0_Nbacmix_10_Nread_Inf_Readlen_50_Npower_0_bacdistflag_0_NoReads.mat']
  
  %auxData.numProcessors = 7;
  auxData.numProcessors = 1;
  auxData.brFlag = 0;
  
  save ~/CS/tmp/run42 auxData
  
  distributeBAC_generalOrSecond('~/CS/tmp/run42');
end


if 1==2
% Compile All POC-related mex files
% mex dna_utils.cpp % this one doesn't work - linker problems
mex calc_all_sites_scores.cpp dna_utils.cpp
mex calc_all_sites_scores_diff_lengths.cpp dna_utils.cpp
mex get_scores_above_threshold.cpp dna_utils.cpp
mex get_scores_above_threshold_diff_lengths.cpp dna_utils.cpp
mex match_closest_vals_c.cpp dna_utils.cpp
mex intervals_intersect_c.cpp dna_utils.cpp hmm_chrom_funcs.cpp
mex extract_sub_kmers.cpp dna_utils.cpp hmm_chrom_funcs.cpp
mex add_noise_to_kmers.cpp dna_utils.cpp hmm_chrom_funcs.cpp

end


if 1==2
  clear
  load ~/tmp/structFor_ffMeanInsteadOfMajority_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_Inf_Noise_0_Correction_0
  load ~/tmp/ffMeanInsteadOfMajority_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_Inf_Noise_0_Correction_0_Nbacmix_10_Nread_Inf_Readlen_50_Npower_05_bacdistflag_0_NoReads

  userdir = getuserdir;
  auxData.saveName = [userdir,'/CS/BAC/ffMeanInsteadOfMajority/ffMeanInsteadOfMajority_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_Inf_Noise_0_Correction_0'];
  auxData.saveName = [userdir,'/CS/BAC/ffMeanInsteadOfMajority/ffMeanInsteadOfMajority_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_Inf_Noise_0_Correction_0'];
  auxData.currDir = [userdir,'/CS/BAC/tmpRuns/ffMeanInsteadOfMajority_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_Inf_Noise_0_Correction_0'];
  auxData.basicSeqNameDir = [userdir,'/CS/BAC/datNoNonACGT/packed64/'];
  auxData.basicSeqKey = [userdir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
  auxData.reads_data_fileName = [userdir,'/tmp/ffMeanInsteadOfMajority_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_Inf_Noise_0_Correction_0_Nbacmix_10_Nread_Inf_Readlen_50_Npower_05_bacdistflag_0_NoReads']
  
  auxData.numProcessors = 7;
  %auxData.numProcessors = 1;
  auxData.brFlag = 0;
  
  save ~/CS/tmp/run1 auxData
  
  distributeBAC_generalOrSecond('~/CS/tmp/run1');

  %ffMeanInsteadOfMajority_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_2_Nreads_Inf_Noise_0_Correction_0
  
end

if 1==2
  clear
  name = '~/CS/tmp/ff1/structFor_ff1_bac_dist_flag_0_Nbac_in_mixture_10_npower_0_readlen_400_numIter_1_Nreads_Inf_Noise_0_Correction_0';
  name1 = find(name=='_');
  name1 = name(name1(1)+1:end);
  load(name)  
  auxData.brFlag = 0;
  auxData.basicSeqNameDir = '~/CS/BAC/datNoNonACGT/packed64/';
  auxData.basicSeqKey = '~/CS/BAC/datNoNonACGT/keyNoNonACGT';
  auxData.reads_data_fileName = ['~/CS/BAC/dataForSim/ff1/',name1,'_Nbacmix_10_Nread_Inf_Readlen_400_Npower_0_bacdistflag_0_NoReads'];
  save(name,'auxData')
  distributeBAC_generalOrSecond(name);
  
  
end