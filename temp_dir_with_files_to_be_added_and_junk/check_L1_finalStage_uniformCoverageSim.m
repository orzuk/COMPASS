function check_L1_finalStage_uniformCoverageSim(basicName,auxDataFileName,resFile)

keyboard
[uniqueReads,uniqueReads_length,auxData_orig] = createReadsFromStructFile(auxDataFileName);

% load L2 results
load(resFile)


%tmpInd = resCell{end}(:,1);%save_store_kp{end};
tmpInd = save_store_kp{end};
[normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOrFourth(auxData_orig.readLength,tmpInd,auxData_orig.basicSeqNameDir,auxData_orig.basicSeqKey);
  
dataIn = struct;
[fracRelevantReads,sumRelevantReads] = currReadsFourth(uniqueReads,uniqueReads_length,values,0,dataIn);

[x2]=runOneGroupOf1000ForCompilationFourth(normalizedBac,fracRelevantReads);    


% solve L1 for the same group
[x1_2_nor] = testL1_2(normalizedBac./max(fracRelevantReads),fracRelevantReads./max(fracRelevantReads));    


% solve L1 for a smaller group
tmpInd_smaller = tmpInd(find(x2>10^-3));
[normalizedBac_smaller values_smaller] = prepareGroupOf1000DistributedSequenceFilesOrFourth(auxData_orig.readLength,tmpInd_smaller,auxData_orig.basicSeqNameDir,auxData_orig.basicSeqKey);
  
dataIn = struct;
[fracRelevantReads_smaller,sumRelevantReads_smaller] = currReadsFourth(uniqueReads,uniqueReads_length,values_smaller,0,dataIn);
[x1_2_nor_smaller] = testL1_2(normalizedBac_smaller./max(fracRelevantReads_smaller),fracRelevantReads_smaller./max(fracRelevantReads_smaller));    

% compare former L2 to current L2
plot(abs(x2(find(x2>10^-5))),abs(resCell{end}(:,2)),'.')


f1 = zeros(auxData_orig.numBACtoConsider,1);
f1(tmpInd_smaller) = x1_2_nor_smaller;

f2 = zeros(auxData_orig.numBACtoConsider,1);
f2(resCell{end}(:,1)) = abs(resCell{end}(:,2));

c = zeros(auxData_orig.numBACtoConsider,1);
c(resCell{1}(:,1)) = abs(resCell{1}(:,2));

mx = max([f1;f2])
subplot(3,1,1)
plot(c,f1,'.')
set(gca,'ylim',[0 mx])
title('L1')
xlabel('c');ylabel('L1')
subplot(3,1,2)
plot(c,f2,'.')
title('L2')
xlabel('c');ylabel('L2')
set(gca,'ylim',[0 mx])
subplot(3,1,3)
plot(f1,f2,'.')
xlabel('L1');ylabel('L2')
set(gca,'ylim',[0 mx])



