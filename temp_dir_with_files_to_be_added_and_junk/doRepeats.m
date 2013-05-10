function [currX,currSumRelevantReads]=doRepeats(uniqueReads,uniqueReads_length,curr_kp,parallelType,basicSaveName,k,userDir,runFileName,auxData)


keyboard
% assumes that curr_kp is ~100000

orig_kp = curr_kp;
tmp_curr_kp = curr_kp;

keepThese = cell(1,auxData.repeatRandomGroups);
for i=1:auxData.repeatRandomGroups
  if length(curr_kp)>auxData.repeatWhenLowerThanThisValue % from 20000 to 150000, do one regular and extract those that are higher than the threshold
    [tmp_currX,tmp_currSumRelevantReads] = iterateParallelDistributedSequenceFilesOr2(uniqueReads,uniqueReads_length,tmp_curr_kp,auxData.basicSeqNameDir,auxData.basicSeqKey,parallelType,[basicSaveName,'_k_',num2str(k+1000)],userDir,[runFileName,'_k_',num2str(k+1000)],auxData);
    
    tmp_keep = find(tmp_currX>auxData.thresholdForCollectingBAC);
    keepThese{i} = tmp_curr_kp(tmp_keep);
    keeThese_X{i} = tmp_currX(tmp_keep);
    tmp_curr_kp(tmp_keep) = [];
  else % lower than ~20000 - do in one batch
    auxData.repeatRandomGroups = auxData.repeatRandomGroups-i+1;
    
    [tmp_currX,tmp_currSumRelevantReads] = calcInBetween(tmp_curr_kp,uniqueReads,uniqueReads_length,parallelType,basicSaveName,k,userDir,runFileName,auxData)
    
    tmp_keep = find(tmp_currX>auxData.thresholdForCollectingBAC);
    keepThese{i} = tmp_curr_kp(tmp_keep);
    keeThese_X{i} = tmp_currX(tmp_keep);
    tmp_curr_kp(tmp_keep) = [];
    break
  end
  
end

if i<auxData.repeatRandomGroups
  keepThese(i+1:auxData.repeatRandomGroups) = [];
  keepThese_X(i+1:auxData.repeatRandomGroups) = [];
end



keep = [];
keep_X = [];
for i=1:length(keepThese)
  if ~isempty(intersect(keep,keepThese{i}))
    disp('problem doRepeats.m');
    pause
  end
  keep = [keep,keepThese{i}];
  keep_X = [keep_X,keeThese_X{i}];
end
% do regular


currX = zeros(1,length(curr_kp));
currX(keep) = keep_X;
currSumRelevantReads = 0;

