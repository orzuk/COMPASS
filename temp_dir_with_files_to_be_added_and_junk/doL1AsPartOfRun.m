function doL1AsPartOfRun(uniqueReads,uniqueReads_length,tmpInd,fileName,auxData)

%keyboard

% rmpath(genpath('/home/csfaculty/shental/Weizmann/ph2users/fenoam/'));addpath('/homes/csfaculty/shental/CS/mFilesBAC/');cd /homes/csfaculty/shental/CS/BAC/cvx; mcc -m doL1AsPartOfRun.m -a prepareGroupOf1000DistributedSequenceFilesOrFourth.m -a currReadsFourth.m -a testL1_2.m -a ./builtins -a ./commands -a ./functions -a ./keywords -a ./lib -a ./sdpt3 -a ./sedumi -a ./sets -a ./structures


[normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOrFourth(auxData.readLength,tmpInd,auxData.basicSeqNameDir,auxData.basicSeqKey);
  
dataIn = struct;
[fracRelevantReads,sumRelevantReads] = currReadsFourth(uniqueReads,uniqueReads_length,values,0,dataIn);

% solve L1 for the same group
[x1_2_nor] = testL1_2(normalizedBac./max(fracRelevantReads),fracRelevantReads./max(fracRelevantReads));    

save([fileName,'_L1Res'],'x1_2_nor','uniqueReads','uniqueReads_length','normalizedBac','values','fracRelevantReads');
