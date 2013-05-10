function [currX,currSumRelevantReads]=repeatedPartitions(uniqueReads,uniqueReads_length,curr_kp,parallelType,basicSaveName,k,userDir,runFileName,auxData)

% split into groups of maximal size of 400
%keyboard

if length(curr_kp)<auxData.repeatWhenLowerThanThisValue
  disp('lower than 20000 - does twice the amount')
  auxData.repeatRandomGroups = 2*auxData.repeatRandomGroups;
end


numRuns = round(length(curr_kp)/auxData.groupSize);

parts = 1:auxData.batchSize:numRuns*auxData.repeatRandomGroups;
if parts(end)~=numRuns*auxData.repeatRandomGroups
  parts(end+1) = numRuns*auxData.repeatRandomGroups+1;
end


clear tmpRepeatRandomGroups;
for i=2:length(parts)
  tmpRepeatRandomGroups(i) = round((parts(i)-parts(i-1))/numRuns);
end
tmpRepeatRandomGroups(1) = [];
if tmpRepeatRandomGroups(end)==0
  tmpRepeatRandomGroups(end) = 1;
end

%keyboard
disp(['total number of repeats: ',num2str(sum(tmpRepeatRandomGroups)),' doRepeats2.m']);

value_forRepeatRandomGroups = zeros(length(curr_kp),sum(tmpRepeatRandomGroups));
ind_col = 0;
for i=1:length(tmpRepeatRandomGroups)
  value_forRepeatRandomGroups(:,ind_col+1:ind_col+tmpRepeatRandomGroups(i)) = calcPartFourth(curr_kp,uniqueReads,uniqueReads_length,parallelType,basicSaveName,k,userDir,runFileName,auxData,tmpRepeatRandomGroups(i));  
  ind_col = ind_col+tmpRepeatRandomGroups(i);
end

mx_majority = mean(value_forRepeatRandomGroups,2);
keepMajority = find(mx_majority>auxData.thresholdForCollectingBAC);

currX = zeros(1,length(curr_kp));
currX(keepMajority) = mx_majority(keepMajority);
currSumRelevantReads = 0;



