function [uniqueReads,uniqueReads_length,auxData] =createReadsFromStructFile(auxDataFileName)




if isdeployed  %load a file which has a variable by the name auxData 
  load(auxDataFileName)
else
  load(auxDataFileName)
end

if ~isfield(auxData,'createSpecificMixtureFlag')
  auxData.createSpecificMixtureFlag = 0;
end


if ~isfield(auxData,'datasetprimers750flag')
  auxData.datasetprimers750flag = 0;
end

if ~isfield(auxData,'createReadsAndQuit')
  auxData.createReadsAndQuit = 0;
end

if ~isfield(auxData,'loadSavedReadsForLargeReadLengthAndManyReads')
  auxData.loadSavedReadsForLargeReadLengthAndManyReads = 0;
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
    disp('ErrorStruct??!! distributeBAC_generalOrFourth.m')
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


%if ~isfield(auxData,'numBACtoConsider')
%  auxData.numBACtoConsider = 455055;
%end

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

if auxData.dataSetPrimers750Flag==0
  if ~isfield(auxData,'basicSeqNameDir') 
    auxData.basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];
  end
  if ~isfield(auxData,'basicSeqKey')
    auxData.basicSeqKey = [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
  end
  disp(['working with database: ',auxData.basicSeqNameDir,'; Is it ok?'])

else
  
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

%keyboard
fn = fieldnames(auxData);
for i=1:length(fn)
    if ~isnumeric(auxData.(fn{i}))
      auxData.(fn{i}) = strrep(auxData.(fn{i}),'/seq/orzuk2/compressed_sensing/metagenomics/next_gen','~');
    end
end    

auxData.userDir = getuserdir;
userDir = getuserdir;
[uniqueReads,uniqueReads_length,auxData] = createReadsData(auxData,userDir);



