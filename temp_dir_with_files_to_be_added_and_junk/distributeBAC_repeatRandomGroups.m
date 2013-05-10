function distributeBAC_repeatRandomGroups(fileToLoad,l)

%%%%%%%%%%%%%%%%%%
load(fileToLoad)
[cX{l},cSumRelevantReads] = iterateParallelDistributedSequenceFilesOr(uniqueReads,uniqueReads_length,curr_kp,auxData.basicSeqNameDir,auxData.basicSeqKey,parallelType,[basicSaveName,'_k_',num2str(k),'_repeat_',num2str(l)],userDir,[runFileName,'_k_',num2str(k),'_repeat_',num2str(l)],auxData);


