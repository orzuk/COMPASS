function extractResRealData_generic(outputName,readDataFileName,resultsDir,auxDataTest)
%keyboard

load([resultsDir,'/',auxDataTest.basicDirNameForFile,'/sol_',outputName,'_noCorrection_',num2str(auxDataTest.readLength)]); 
load([readDataFileName],'uniqueReads','uniqueReads_length')

find(abs(found{end})>10^-3);

userDir = getuserdir;

if strcmp(auxDataTest.datasetName,'fullDatabse')
  basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];
  basicSeqKey= [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
else
  basicSeqNameDir = [userDir,'/CS/BAC/',auxDataTest.datasetName,'/datNoNonACGT/packed64/'];
  basicSeqKey= [userDir,'/CS/BAC/',auxDataTest.datasetName,'/datNoNonACGT/keyNoNonACGT_',auxDataTest.datasetName];
  
end

tmpInd = find(abs(found{end})>10^-3);
[normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOrFourth(auxDataTest.readLength,tmpInd,basicSeqNameDir,basicSeqKey);

dataIn = struct;
[fracRelevantReads,sumRelevantReads] = currReadsFourth(uniqueReads,uniqueReads_length,values,0,dataIn);
%keyboard 
dep = found{end}(tmpInd);

[x1_2_nor] = testL1_2(normalizedBac./max(fracRelevantReads),fracRelevantReads./max(fracRelevantReads));    
%[x1_2_no_nor] = testL1_2(normalizedBac,fracRelevantReads);    


numBACtoConsider = auxDataTest.numBACtoConsider;
f_l1 = zeros(numBACtoConsider,1);
f_l1(tmpInd) = x1_2_nor;

f_l2 = zeros(numBACtoConsider,1);
f_l2(tmpInd) = abs(found{end}(tmpInd));

figure(1);clf
if 1==2
clf
plot(c,f_l1,'b.')
hold on
plot(c,f_l2,'kx')
legend('l1','l2')
title(inputName)
%print('-dpdf',['~/CS/BAC/',inputName])  
figure(2);clf
plot(c,f_l1,'kx')
line([0,max(c)],[0 max(c)])
xlabel('correct');ylabel('l1')
figure(3);clf
plot(c,f_l2,'kx')
line([0,max(c)],[0 max(c)])
xlabel('correct');ylabel('l2')
end

figure(4);clf
plot(f_l1,f_l2,'.')
%line([0,max(c)],[0 max(c)])
xlabel('l1');ylabel('l2')
%title(fileName)



%keyboard
save([resultsDir,'/',auxDataTest.basicDirNameForFile,'/res_',outputName],'f_l1','f_l2','normalizedBac','values','x1_2_nor','tmpInd','fracRelevantReads')
