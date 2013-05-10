function [currX,currSumRelevantReads]=calcInBetweenFourth(curr_kp,uniqueReads,uniqueReads_length,parallelType,basicSaveName,k,userDir,runFileName,auxData)

disp(['running with: ',num2str(length(curr_kp)),'; repeating runs: ',num2str(auxData.repeatRandomGroups),' calcInBetween.m'])
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

auxData.keepOriginalOrderFlag = 1;
auxData.forcePartFromOutsideFlag = 1;
auxData.partFromOutside = currPart;
%keyboard
[cX,cSumRelevantReads] = iterateParallelDistributedSequenceFilesOrFourth(uniqueReads,uniqueReads_length,curr_kp_forRepeatRandomGroups,auxData.basicSeqNameDir,auxData.basicSeqKey,parallelType,[basicSaveName,'_k_',num2str(k+2000)],userDir,[runFileName,'_k_',num2str(k+2000)],auxData);

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
  disp('average the repeatRandomGroups in distributeBAC_generalOrThird.m. Does it work well?')
  mx_majority = mean(value_forRepeatRandomGroups,2);
  keepMajority = find(mx_majority>auxData.thresholdForCollectingBAC);
end

currX = zeros(1,length(curr_kp));
currX(keepMajority) = mx_majority(keepMajority);
currSumRelevantReads = 0;
