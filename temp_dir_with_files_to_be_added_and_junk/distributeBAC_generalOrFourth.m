function distributeBAC_generalOrFourth(auxDataFileName)

% the same as distributeBAC_general but also does infinite reads. allows differnet normalizations for the number of reads

% reads are red
% can load specific reads - just like needed in similar BAC! just make sure the reads are called red

% create stand alone run files and then run them

% do queue L1
%rmpath(genpath('/home/csfaculty/shental/Weizmann/ph2users/fenoam/'));rmpath('/home/csfaculty/shental/matlab');cd ~/CS/mFilesBAC;addpath('~/CS/mFilesBAC');  mcc -v -m distributeBAC_generalOrFourth  -d /homes/csfaculty/shental/CS/BAC/cvx -a grep.m  -a currReadsFourth.m -a findParametersReads.m -a fun_prep_data_for_simulations_loadSpecificBAC2OrFourth.m -a get_duplicates.m -a get_duplicates2.m -a getuserdir.m -a iterateParallelDistributedSequenceFilesOrFourth.m  -a myHistForBAC.m -a myUniqueUINT8.m -a randi.m -a newReadsBasedSeed2Or.m   -a prepareReadsForRandomSequences.m -a solveForGroupOrFourth.m   -a prepareGroupOf1000DistributedSequenceFilesOrFourth  -a runOneGroupOf1000ForCompilationFourth.m -a set_solution_bounds.m -a doRepeatsFourth.m -a calcPartFourth.m -a createReadsData -a CreateMothurDistNoam2.m -a doL1AsPartOfRun.m -a CreateMothurDistNoamL1L2.m


% do local L1
% compile_distributeBAC_generalOrFourth



if isdeployed  %load a file which has a variable by the name auxData 
  load(auxDataFileName)
else
  load(auxDataFileName)
end
%keyboard

% change 
if ~isempty(findstr(getenv('HOSTNAME'),'openu'))
  auxData.reads_data_fileName = strrep(auxData.reads_data_fileName,'/seq/orzuk2/compressed_sensing/metagenomics/next_gen',getuserdir); 
  auxData.userDir = getuserdir;
  auxData.basicSeqNameDir = strrep(auxData.basicSeqNameDir,'/seq/orzuk2/compressed_sensing/metagenomics/next_gen',getuserdir); 
  auxData.basicSeqKey = strrep(auxData.basicSeqKey,'/seq/orzuk2/compressed_sensing/metagenomics/next_gen',getuserdir); 
  auxData.brFlag = 0;
  %keyboard
end

if ~isfield(auxData,'doMothurFlag') % Solution evaluation using hte Mothur flag
  auxData.doMothurFlag = 0;
end


if ~isfield(auxData,'doL1Flag')  % Run robust L1 optimization for real data 
  auxData.doL1Flag = 0;
end


if ~isfield(auxData,'createSpecificMixtureFlag') % Set identities of species in mixture 
  auxData.createSpecificMixtureFlag = 0;
end


if ~isfield(auxData,'datasetprimers750flag') % Run on 750 primers (default 'NO': full 16S)
  auxData.datasetprimers750flag = 0;
end

if ~isfield(auxData,'createReadsAndQuit')  % Only simulation 
  auxData.createReadsAndQuit = 0;
end

if ~isfield(auxData,'loadSavedReadsForLargeReadLengthAndManyReads')  % load reads from file ??? 
  auxData.loadSavedReadsForLargeReadLengthAndManyReads = 0;
end


if ~isfield(auxData,'correctMeasurementForReadErrorsFlag')  % correct read errors
  auxData.correctMeasurementForReadErrorsFlag = 0;
end

if ~isfield(auxData,'batchSize') % number of runs per run in cluster - timewise!
  auxData.batchSize = 400;
end

if ~isfield(auxData,'upperLimitForRepeat')  % Algorithmic parameter: perform repetition below this value 
  auxData.upperLimitForRepeat = 150000;
end

if ~isfield(auxData,'repeatWhenLowerThanThisValue') % Algorithmic parameter: perform repetition above (??) this value
  auxData.repeatWhenLowerThanThisValue = 20000;
end

if ~isfield(auxData,'saveFullDataFlag') % Save format (always .mat file) 
  auxData.saveFullDataFlag = 0;
end

if ~isfield(auxData,'realReadsFromFileFlag')  % Real data (we perform additional pre-processing step) 
  auxData.realReadsFromFileFlag = 0;
end

if ~isfield(auxData,'addNoiseFlag') % Add noise to simulated reads 
  auxData.addNoiseFlag = 0;
end

if auxData.addNoiseFlag  % load read-error parameters
  if ~isfield(auxData,'ErrorStruct')
    disp('ErrorStruct??!! distributeBAC_generalOrFourth.m')
    return
  end
else
  auxData.ErrorStruct = [];
end

if ~isfield(auxData,'forcePartFromOutsideFlag') % Algorithmic parameter: set how to divide bacteria to blocks
  auxData.forcePartFromOutsideFlag = 0;
end

if ~isfield(auxData,'preprocessByExistanceFlag') % ??????
  auxData.preprocessByExistanceFlag = 0;
end
if ~isfield(auxData,'repeatRandomGroups') % Algorithmic parameter: how many repeat random selections of groups performed 
  auxData.repeatRandomGroups = 1;
end

if ~isfield(auxData,'inifiniteNumberOfReadsFlag') % Simulate infinite reads (no noise). 
  auxData.inifiniteNumberOfReadsFlag = 0;
end

if ~isfield(auxData,'keepOriginalOrderFlag') % keepOriginalOrderFlag for checking in the original order
  auxData.keepOriginalOrderFlag = 0;
end


if ~isfield(auxData,'randomizeSeedAfterCreateReadsFlag') % for Matlab random number generator: randomizeSeedAfterCreateReadsFlag for considering a certain bacteria from a group
  auxData.randomizeSeedAfterCreateReadsFlag = 0;
end



% add to auxData stuff
if ~isfield(auxData,'groupSize') % Algorithmic parameter: partition into groups of groupSize bacteria (this is the block size) 
  auxData.groupSize = 1000;
end

if ~isfield(auxData,'smallestSetOfCollected') % Algorithmic parameter: Stopping criteria for divide-and-conquor. After collecting the bacteria from each group - continue iterating id their number is smaller than smallestSetOfCollected
  auxData.smallestSetOfCollected = 1000;
end

if ~isfield(auxData,'thresholdForCollectingBAC') % Algorithmic parameter: Frequency threshold for keeping a certain bacteria from a group. (Threshold to zero all bacteria below this value) 
  auxData.thresholdForCollectingBAC = 1e-3;
end


%if ~isfield(auxData,'numBACtoConsider')
%  auxData.numBACtoConsider = 455055;
%end

if ~isfield(auxData,'createReadsFlag')  % again: only create reads ??? 
  auxData.createReadsFlag = 0;
end

%%%%%% %end add to auxData stuff


if ~isfield(auxData,'brFlag') % Run in broad cluster 
  auxData.brFlag = 0;
end

if ~isfield(auxData,'queueName') % queue to run jobs 
  auxData.queueName = [];
end

if ~isfield(auxData,'firstFlag') % ?????? 
  auxData.firstFlag = 0;
end

if auxData.brFlag  % run on broad cluster
    userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
    disp('run in Br')
    parallelType = 'cluster';
else  % run on different computer: only one computer 
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

if auxData.dataSetPrimers750Flag==0  % 750 primer flag 
  if ~isfield(auxData,'basicSeqNameDir') 
    auxData.basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];
  end
  if ~isfield(auxData,'basicSeqKey')
    auxData.basicSeqKey = [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
  end
  disp(['working with database: ',auxData.basicSeqNameDir,'; Is it ok?'])

else  % Full 16S flag 
  
  if ~isfield(auxData,'basicSeqNameDir') 
    auxData.basicSeqNameDir = [userDir,'/CS/BAC/primers750/datNoNonACGT/packed64/'];
  end
  if ~isfield(auxData,'basicSeqKey')
    auxData.basicSeqKey = [userDir,'/CS/BAC/primers750/datNoNonACGT/keyNoNonACGT_primers750'];
  end
  disp(['working with database: ',auxData.basicSeqNameDir,'; Is it ok?'])

  
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end load variables

% Write parameters to screen 
fn = fieldnames(auxData);
fprintf('%%%%%%%%%%%\n');
fprintf('data:\n');
for i=1:length(fn)
  if ~strcmp(fn{i},'ErrorStruct') & ~strcmp(fn{i},'readsForCorrectReads')
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

if ~isdeployed && auxData.firstFlag==1   % Add additional packages to matlab path 
  addpath([userDir,'/CS/mFilesBAC/']);
  addpath(genpath([userDir,'/CS/BAC/cvx/']))
end


%keyboard

%unix(['mkdir ',auxData.currDir]);
runFileName = [auxData.currDir,'/run_',auxData.tmpFileName]; % Prepare jobs 
basicSaveName = [auxData.currDir,'/',auxData.tmpFileName,'_dat'];
%keyboard

% create the reads or read them from file
if auxData.realReadsFromFileFlag==0
  disp('using simulated reads')
  if auxData.inifiniteNumberOfReadsFlag==0
    disp('finite number of reads')
    
    if auxData.loadSavedReadsForLargeReadLengthAndManyReads==0
      disp('run standard. create reads. distributeBAC_generalOrFourth.m')
      % change 4.1.12 - output the auxData, Necessary for the myCorrectReadsNew
      [uniqueReads,uniqueReads_length,auxData] = createReadsData(auxData,userDir);
      
    else % auxData.loadSavedReadsForLargeReadLengthAndManyReads=1
      
      
      if auxData.createReadsAndQuit==1
        disp('creating reads and quit. distributeBAC_generalOrFourth.m')
        % change 16.1.12
        %[uniqueReads,uniqueReads_length] = createReadsData(auxData,userDir);
        [uniqueReads,uniqueReads_length,auxData] = createReadsData(auxData,userDir);
        % end change 16.1.12
        
        auxData.createReadsAndQuit=0;
        auxData.loadSavedReadsForLargeReadLengthAndManyReads=1;
        %keyboard
        
        createReadsIndex = findstr(auxDataFileName,'_createReads');
        auxDataFileName(createReadsIndex:end) = [];
        
         if ~isempty(findstr(getenv('HOSTNAME'),'openu'))
           auxData.reads_data_fileName = strrep(auxData.reads_data_fileName,getuserdir,'/seq/orzuk2/compressed_sensing/metagenomics/next_gen'); 
           auxData.userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
           auxData.basicSeqNameDir = strrep(auxData.basicSeqNameDir,getuserdir,'/seq/orzuk2/compressed_sensing/metagenomics/next_gen'); 
           auxData.basicSeqKey = strrep(auxData.basicSeqKey,getuserdir,'/seq/orzuk2/compressed_sensing/metagenomics/next_gen'); 
           auxData.brFlag = 1;
           %keyboard
         end 
        
        
         save(auxDataFileName,'auxData','uniqueReads','uniqueReads_length');  
        return
      else
        disp('read should already be loaded. createReadsAndQuit=1. distributeBAC_generalOrFourth.m')
        
      end
      
      
    end
      
    %%%%%%%%%%%%%%%%%
    
  else
    disp('using an infinite number of reads');
    
    load(auxData.reads_data_fileName,'correctWeight','ind_bac_in_mix')
    
    % for infinite 
    
    [uniqueReads,uniqueReads_length,auxData.fracRelevantReadsForInfinity] = createReadsForInfiniteNumberOrFourth(ind_bac_in_mix,correctWeight,auxData.readLength,auxData.basicSeqNameDir,auxData.basicSeqKey);
    
    clear ind_bac_in_mix correctWeight
  end

else % real reads
  disp('using reads from file')
  load(auxData.reads_data_fileName) % uniqueReads uniqueReads_length
  
  
  
  
end % load reads

%keyboard

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
      disp('running one random split. iterateParallelDistributedSequenceFilesOrFourth.m ')
      [currX,currSumRelevantReads] = iterateParallelDistributedSequenceFilesOrFourth(uniqueReads,uniqueReads_length,curr_kp,auxData.basicSeqNameDir,auxData.basicSeqKey,parallelType,[basicSaveName,'_k_',num2str(k)],userDir,[runFileName,'_k_',num2str(k)],auxData);
    else
      % less than 150000. when lower than 50000 does in one batch
      [currX,currSumRelevantReads] = doRepeatsFourth(uniqueReads,uniqueReads_length,curr_kp,parallelType,basicSaveName,k,userDir,runFileName,auxData);
      %what about the value of k?
    end
    
    
  else % regular
    disp('running one random split. iterateParallelDistributedSequenceFilesOrFourth.m ')
    [currX,currSumRelevantReads] = iterateParallelDistributedSequenceFilesOrFourth(uniqueReads,uniqueReads_length,curr_kp,auxData.basicSeqNameDir,auxData.basicSeqKey,parallelType,[basicSaveName,'_k_',num2str(k)],userDir,[runFileName,'_k_',num2str(k)],auxData);
  
  end
  
  
  
  
  next_kp = curr_kp(find(currX>auxData.thresholdForCollectingBAC));
  store_kp{k} = curr_kp;
  store_X{k} = currX;
  store_sumRelevantReads{k} = currSumRelevantReads;
  k = k +1;
  
  curr_kp = next_kp;
  
  if length(store_kp{k-1})==length(curr_kp)
    disp('did not reduce size, so breaks. distributeBAC_general3.m')
    save([auxData.saveName,'_break'],'store_kp','store_X','store_sumRelevantReads')
    break
  end
end
if length(store_kp{k-1})~=length(curr_kp) % size was reduced
  store_kp{k} = curr_kp;
  store_X{k} = [];
end

%keyboard
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

% change 12.12.11
if isempty(do) && lengths_kp(end)<3000
  do = length(lengths_kp);
end
% end change 12.12.11
%keyboard
% deal with the case that n-1 is larger than auxData.smallestSetOfCollected*2 and the n'th is really small

if auxData.realReadsFromFileFlag==0
  if isempty(find(lengths_kp(do)>100))
    do = do(1)-1;
    disp(['running with larger group than ',num2str(auxData.smallestSetOfCollected*2),' of size: ',num2str(lengths_kp(do)),' distributeBAC_general3.m'])
  end
else % real reads
   %real reads - take all do
  
end

% solve - regular but the negative solution are the ones which are unique_inds
auxData.solveOr = 1;
auxData.do_lp_flag = 0;
auxData.numProcessors = 1;
auxData.moveDependentToNextStageFlag = 0;
%auxData.keepOriginalOrderFlag==1;
for i=do
  curr_kp = store_kp{i};
  auxData.groupSize = length(curr_kp)-1;
  [currX,currSumRelevantReads] = iterateParallelDistributedSequenceFilesOrFourth(uniqueReads,uniqueReads_length,curr_kp,auxData.basicSeqNameDir,auxData.basicSeqKey,parallelType,[basicSaveName,'_final_',num2str(i)],userDir,[runFileName,'_final_',num2str(i)],auxData);
  found{i} = zeros(1,auxData.numBACtoConsider);
  found{i}(curr_kp) = currX;
  found_sumRelevantReads(i) = currSumRelevantReads;
end

if auxData.realReadsFromFileFlag==0
  if auxData.saveFullDataFlag
    load(auxData.reads_data_fileName,'correctWeight','ind_bac_in_mix')  
    disp('standard way of saving')
    save(auxData.saveName,'found','found_sumRelevantReads','store_kp','store_X','store_sumRelevantReads','correctWeight','ind_bac_in_mix')
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
    %%%%%%%
    % save mothur - in the same directory
    
%    if auxData.doL1Flag==0
%      CreateMothurDistNoam2(auxData.saveName,0,auxData); % does mothur on 10^-5 of the resCell{2}
%    end
      
    
  end
else % real reads
  save(auxData.saveName,'found','found_sumRelevantReads','store_kp','store_X','store_sumRelevantReads')
end
%keyboard
% save L1 data
if auxData.doL1Flag
  disp('notice: the matrix without weights!!!!!')

  tmpInd = store_kp{end};
  fileNameL1Res = [auxData.saveName,'_L1Res'];
  basicAllocationFileName = [auxData.saveName,'_L1data'];
  save(basicAllocationFileName,'uniqueReads','uniqueReads_length','tmpInd','auxData','fileNameL1Res')
  
  % run in queue
  if 1==2
    %runFileName = [auxData.saveName,'_tmpRunL1'];
    %doL1AsPartOf_distribute(basicAllocationFileName,runFileName,fileNameL1Res,auxData)
  end
  
  % run here
  x1_2_nor = addL1toL2_res(uniqueReads,uniqueReads_length,tmpInd,fileNameL1Res,auxData);
  
  % do mothur of both L1 and L2 - mothur is done on resCell{end} which is the L1 result
%  if auxData.realReadsFromFileFlag==0 %simulated reads
%    CreateMothurDistNoamL1L2(auxData.saveName,tmpInd',x1_2_nor,resCell{end}(:,1),resCell{end}(:,2),correctWeight,auxData);
     
%    README = 'save_store_kp{end} is the basis of both results: L2 takes resCell{end}(:,2)> 10^-5 and L1 calculates L1 over save_store_kp{end}'
%  end
README = 'save_store_kp{end} is the basis of both results: L2 takes resCell{end}(:,2)> 10^-5 and L1 calculates L1 over save_store_kp{end}';  
  
  if isempty(find(isnan(x1_2_nor))) % was solved
   % once sure it's ok - del the file
   unix(['rm -f ',basicAllocationFileName,'.mat']);
  end
  
    
end



% delete tmpRuns
if exist('found')
  %unix(['rm -f ',auxData.currDir,'/',auxData.tmpFileName,'_*']);
  unix(['rm -fr ',auxData.currDir,'/']);
  %unix(['rm -f ',auxData.currDir,'/run_',auxData.tmpFileName,'_*']);
end


if 1==2
  clear
  load /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/canonicalFormWithL1/canonicalFormWithL1_bac_dist_flag_0_Nbac_in_mixture_1000_npower_5_readlen_100_numIter_6_Nreads_1000000_Noise_0_Correction_0_L1data
  fileNameL1Res = [auxData.saveName,'_L1Res'];
  save /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/canonicalFormWithL1/canonicalFormWithL1_bac_dist_flag_0_Nbac_in_mixture_1000_npower_5_readlen_100_numIter_6_Nreads_1000000_Noise_0_Correction_0_L1data

  clear
  load /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/canonicalFormWithL1/canonicalFormWithL1_bac_dist_flag_0_Nbac_in_mixture_1000_npower_5_readlen_100_numIter_6_Nreads_1000000_Noise_0_Correction_0_L1data
  basicAllocationFileName = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/canonicalFormWithL1/canonicalFormWithL1_bac_dist_flag_0_Nbac_in_mixture_1000_npower_5_readlen_100_numIter_6_Nreads_1000000_Noise_0_Correction_0_L1data';
  runFileName = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/tmpRuns/canonicalFormWithL1/canonicalFormWithL1_bac_dist_flag_0_Nbac_in_mixture_1000_npower_5_readlen_100_numIter_6_Nreads_1000000_Noise_0_Correction_0/junkt';
  fileNameL1Res = [auxData.saveName,'_L1Res'];
  doL1AsPartOf_distribute(basicAllocationFileName,runFileName,fileNameL1Res,auxData)

end


if 1==2
  clear
  load /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/canonicalFormWithL1/canonicalFormWithL1_bac_dist_flag_0_Nbac_in_mixture_1000_npower_5_readlen_100_numIter_6_Nreads_1000000_Noise_0_Correction_0_L1data
  
  
end


