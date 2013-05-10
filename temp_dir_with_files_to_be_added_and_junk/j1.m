tmpInd = reshape(curr_kp_forRepeatRandomGroups,7950,10);

d{1} = curr_kp_forRepeatRandomGroups(1:7950);
%d{1} = tmpInd(:,1);
find(d{1}==71281)

%e{1} = d{1}(6001:7950);
e{1} = d{1}(6001:7000);
find(e{1}==71281)
[x,sumRelevantReads]=solveForGroupOr(1,1,auxData,[],[],e,auxData.basicSeqNameDir,auxData.basicSeqKey,auxData.readLength,uniqueReads,uniqueReads_length);

e{1} = d{1}(6001:7950);
find(e{1}==71281)
[y,sumRelevantReads]=solveForGroupOr(1,1,auxData,[],[],e,auxData.basicSeqNameDir,auxData.basicSeqKey,auxData.readLength,uniqueReads,uniqueReads_length);


find(curr_kp==71281)

