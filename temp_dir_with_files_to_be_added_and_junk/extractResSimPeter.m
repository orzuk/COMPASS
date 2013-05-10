function extractResSimPeter(outputName,inputName,fileName,readLength)
%keyboard

load(['~/CS/BAC/Peter/sol_',outputName,'_',outputName,'_noCorrection_',num2str(readLength),'.mat']); 
%load(['~/CS/BAC/Peter/',inputName],'freq','bact')
load(['~/CS/BAC/Peter/',fileName])


find(abs(found{end})>10^-3)

userDir = getuserdir;
basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];
basicSeqKey = [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];

tmpInd = find(abs(found{end})>10^-3);
[normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOrFourth(readLength,tmpInd,basicSeqNameDir,basicSeqKey);
  
  
dataIn = struct;
[fracRelevantReads,sumRelevantReads] = currReadsFourth(uniqueReads,uniqueReads_length,values,0,dataIn);

[x1_2_nor] = testL1_2(normalizedBac./max(fracRelevantReads),fracRelevantReads./max(fracRelevantReads));    
[x1_2_no_nor] = testL1_2(normalizedBac,fracRelevantReads);    

%keyboard
f_l1 = zeros(410849,1);
f_l1(tmpInd) = x1_2_nor;

f_l2 = zeros(410849,1);
f_l2(tmpInd) = abs(found{end}(tmpInd));


plot(f_l1,f_l2,'.')
%line([0,max(c)],[0 max(c)])
xlabel('l1');ylabel('l2')
title(fileName)
%keyboard
save(['~/CS/BAC/Peter/simRes/res_',inputName],'f_l1','f_l2','normalizedBac','values','x1_2_nor','tmpInd','fracRelevantReads')


