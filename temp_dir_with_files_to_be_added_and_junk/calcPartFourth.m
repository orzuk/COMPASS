function value_forRepeatRandomGroups=calcPartFourth(curr_kp,uniqueReads,uniqueReads_length,parallelType,basicSaveName,k,userDir,runFileName,auxData,tmpRepeatRandomGroups)

disp(['running with: ',num2str(length(curr_kp)),'; repeating runs: ',num2str(tmpRepeatRandomGroups),' calcInBetween.m'])
% assumes curr_kp is sorted

partTemplate = 1:auxData.groupSize:length(curr_kp);
partTemplate(end) = length(curr_kp)+1;

curr_perm = randperm(length(curr_kp));
curr_kp_forRepeatRandomGroups = curr_kp(curr_perm);
currPart = partTemplate;

%keyboard
z = zeros(length(curr_kp),tmpRepeatRandomGroups);
z((curr_perm),1) = 1:length(curr_kp);

for l=2:tmpRepeatRandomGroups
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

if tmpRepeatRandomGroups==1
  value_forRepeatRandomGroups = value_forRepeatRandomGroups';
end

