function runGeneralWithLoadBalance(auxData)

% last change 13.1 - add CorrectReadsNew to crGeneralFourth

% change 14.1 - increase numCPU
if ~isfield(auxData,'numProcessors')
  auxData.numProcessors = 10;
end

defValues = struct;
defValues.bdf = 0;
defValues.groupSize = 1000;
defValues.thresholdForCollectingBAC = 10^-3;

auxData = checkValues(auxData,defValues);

printStructValues(auxData,'auxData in runGeneral');

%keyboard


userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
addpath([userDir,'/CS/mFilesBAC'])
unix(['mkdir ',userDir,'/CS/tmp/',auxData.outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/tmpRuns/',auxData.outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/',auxData.outFileDirName])
prefix = auxData.outFileDirName;
brFlag = 1;
% using load balance
crGeneralFourthWithLoadBalance(prefix,auxData.outFileDirName,auxData.numProcessors,auxData.bdf,auxData.nb,auxData.pwr,auxData.numIter,auxData.nr,auxData.addNoiseFlagInput,auxData.correctMeasurementForReadErrorsFlagInput,auxData.readLength,'hour','week',150000,50000,400,brFlag,1,auxData.thresholdForCollectingBAC,auxData.groupSize,auxData.dataSetPrimers750Flag,auxData.CorrectReadsNew,auxData.doL1Flag,auxData.fileNameForFindNumber);
crGeneralFourthWithLoadBalance(prefix,auxData.outFileDirName,auxData.numProcessors,auxData.bdf,auxData.nb,auxData.pwr,auxData.numIter,auxData.nr,auxData.addNoiseFlagInput,auxData.correctMeasurementForReadErrorsFlagInput,auxData.readLength,'hour','week',150000,50000,400,brFlag,0,auxData.thresholdForCollectingBAC,auxData.groupSize,auxData.dataSetPrimers750Flag,auxData.CorrectReadsNew,auxData.doL1Flag,auxData.fileNameForFindNumber);


listName_noCorrection = dir([userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*_noCorrection*']);
%keyboard
name = ['tr1_2_',auxData.outFileDirName,'_noCorrection.txt']
fid = fopen([userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/',name],'w');
fprintf(fid,'#!/bin/bash \n');

k=1;
for i=1:length(listName_noCorrection)
  w = [listName_noCorrection(i).name,' ;sleep 20;'];
  fprintf(fid,'./%s\n',w);
 
end
fclose(fid);

unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/',name])

%%%%%%%%
% with correction
% sum the number of CPU

listName_withCorrection = dir([userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*_withCorrection*']);
%keyboard
name = ['tr1_2_',auxData.outFileDirName,'_withCorrection.txt']
fid = fopen([userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/',name],'w');
fprintf(fid,'#!/bin/bash \n');

k=1;
for i=1:length(listName_withCorrection)
  w = [listName_withCorrection(i).name,' ;sleep 20;'];
  fprintf(fid,'./%s\n',w);
 
end
fclose(fid);

unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/',name])


