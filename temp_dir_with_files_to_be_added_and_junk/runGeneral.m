function runGeneral(auxData)

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
crGeneralFourth(prefix,auxData.outFileDirName,auxData.numProcessors,auxData.bdf,auxData.nb,auxData.pwr,auxData.numIter,auxData.nr,auxData.addNoiseFlagInput,auxData.correctMeasurementForReadErrorsFlagInput,auxData.readLength,'hour','week',150000,50000,400,brFlag,1,auxData.thresholdForCollectingBAC,auxData.groupSize,auxData.dataSetPrimers750Flag,auxData.CorrectReadsNew);
crGeneralFourth(prefix,auxData.outFileDirName,auxData.numProcessors,auxData.bdf,auxData.nb,auxData.pwr,auxData.numIter,auxData.nr,auxData.addNoiseFlagInput,auxData.correctMeasurementForReadErrorsFlagInput,auxData.readLength,'hour','week',150000,50000,400,brFlag,0,auxData.thresholdForCollectingBAC,auxData.groupSize,auxData.dataSetPrimers750Flag,auxData.CorrectReadsNew);

[junk,list36] = unix(['grep -l "readlen_36" ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*']);
[junk,list25] = unix(['grep -l "readlen_25" ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*']);
[junk,list10] = unix(['grep -l "readlen_10" ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*']);
[junk,list50] = unix(['grep -l "readlen_50" ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*']);
[junk,list100] = unix(['grep -l "readlen_100" ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*']);
[junk,list400] = unix(['grep -l "readlen_400" ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*']);

l10 = extract_listNames(list10);
l25 = extract_listNames(list25);
l36 = extract_listNames(list36);
l50 = extract_listNames(list50);
l100 = extract_listNames(list100);
l400 = extract_listNames(list400);

l = [l400';l100';l50';l36';l25';l10'];
%pos = [60*ones(length(l400),1);90*ones(length(l100),1);4*60*ones(length(l50),1)];
pos = [120*ones(length(l400),1);180*ones(length(l100),1);4*60*ones(length(l50),1);4*60*ones(length(l36),1);4*60*ones(length(l25),1);;4*60*ones(length(l10),1);];

dup = [];
for i=1:length(l)-1
  for j=i+1:length(l)
    if length(l{i})==length(l{j})
      if strcmp(l{i},l{j})
        dup = [dup,j];
      end
    end
  end
end
l(dup) = []; 
pos(dup) = [];


name = 'tr1_2.txt'
fid = fopen([userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/',name],'w');
fprintf(fid,'#!/bin/bash \n');

k=1;
for i=1:length(l)
  w = [l{i},' ;sleep 20;'];
  fprintf(fid,'%s\n',w);
  if mod(k,5)==0
    fprintf(fid,'sleep %s;\n',[num2str(pos(i)),'m']);
  end
  k = k+1;
end
fclose(fid);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/',name])



