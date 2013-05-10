function runGeneralSplitCorrection(auxData)

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
crGeneralFourth_splitCorrection(prefix,auxData.outFileDirName,auxData.numProcessors,auxData.bdf,auxData.nb,auxData.pwr,auxData.numIter,auxData.nr,auxData.addNoiseFlagInput,auxData.correctMeasurementForReadErrorsFlagInput,auxData.readLength,'hour','week',150000,50000,400,brFlag,1,auxData.thresholdForCollectingBAC,auxData.groupSize,auxData.dataSetPrimers750Flag,auxData.CorrectReadsNew,auxData.doL1Flag);
crGeneralFourth_splitCorrection(prefix,auxData.outFileDirName,auxData.numProcessors,auxData.bdf,auxData.nb,auxData.pwr,auxData.numIter,auxData.nr,auxData.addNoiseFlagInput,auxData.correctMeasurementForReadErrorsFlagInput,auxData.readLength,'hour','week',150000,50000,400,brFlag,0,auxData.thresholdForCollectingBAC,auxData.groupSize,auxData.dataSetPrimers750Flag,auxData.CorrectReadsNew,auxData.doL1Flag);

[junk,list35_noCorrection] = unix(['grep -l "readlen_35_" ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*_noCorrection*']);
[junk,list75_noCorrection] = unix(['grep -l "readlen_75_" ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*_noCorrection*']);

[junk,list36_noCorrection] = unix(['grep -l "readlen_36_" ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*_noCorrection*']);
[junk,list25_noCorrection] = unix(['grep -l "readlen_25_" ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*_noCorrection*']);
[junk,list10_noCorrection] = unix(['grep -l "readlen_10_" ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*_noCorrection*']);
[junk,list50_noCorrection] = unix(['grep -l "readlen_50_" ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*_noCorrection*']);
[junk,list100_noCorrection] = unix(['grep -l "readlen_100_" ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*_noCorrection*']);
[junk,list150_noCorrection] = unix(['grep -l "readlen_150_" ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*_noCorrection*']);
[junk,list200_noCorrection] = unix(['grep -l "readlen_200_" ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*_noCorrection*']);

[junk,list400_noCorrection] = unix(['grep -l "readlen_400_" ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*_noCorrection*']);

[junk,list35_withCorrection] = unix(['grep -l "readlen_35_" ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*_withCorrection*']);
[junk,list75_withCorrection] = unix(['grep -l "readlen_75_" ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*_withCorrection*']);
[junk,list36_withCorrection] = unix(['grep -l "readlen_36_" ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*_withCorrection*']);
[junk,list25_withCorrection] = unix(['grep -l "readlen_25_" ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*_withCorrection*']);
[junk,list10_withCorrection] = unix(['grep -l "readlen_10_" ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*_withCorrection*']);
[junk,list50_withCorrection] = unix(['grep -l "readlen_50_" ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*_withCorrection*']);
[junk,list100_withCorrection] = unix(['grep -l "readlen_100_" ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*_withCorrection*']);
[junk,list150_withCorrection] = unix(['grep -l "readlen_150_" ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*_withCorrection*']);
[junk,list200_withCorrection] = unix(['grep -l "readlen_200_" ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*_withCorrection*']);

[junk,list400_withCorrection] = unix(['grep -l "readlen_400_" ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/list*_withCorrection*']);


l10_noCorrection = extract_listNames(list10_noCorrection);
l25_noCorrection = extract_listNames(list25_noCorrection);
l35_noCorrection = extract_listNames(list35_noCorrection);
l36_noCorrection = extract_listNames(list36_noCorrection);
l50_noCorrection = extract_listNames(list50_noCorrection);
l75_noCorrection = extract_listNames(list75_noCorrection);
l100_noCorrection = extract_listNames(list100_noCorrection);
l150_noCorrection = extract_listNames(list150_noCorrection);
l200_noCorrection = extract_listNames(list200_noCorrection);
l400_noCorrection = extract_listNames(list400_noCorrection);

l10_withCorrection = extract_listNames(list10_withCorrection);
l25_withCorrection = extract_listNames(list25_withCorrection);
l35_withCorrection = extract_listNames(list35_withCorrection);
l36_withCorrection = extract_listNames(list36_withCorrection);
l50_withCorrection = extract_listNames(list50_withCorrection);
l75_withCorrection = extract_listNames(list75_withCorrection);
l100_withCorrection = extract_listNames(list100_withCorrection);
l150_withCorrection = extract_listNames(list150_withCorrection);
l200_withCorrection = extract_listNames(list200_withCorrection);
l400_withCorrection = extract_listNames(list400_withCorrection);


l_noCorrection = [l35_noCorrection';l75_noCorrection';l150_noCorrection';l400_noCorrection';l100_noCorrection';l200_noCorrection';l50_noCorrection';l36_noCorrection';l25_noCorrection';l10_noCorrection'];
l_withCorrection = [l35_withCorrection';l75_withCorrection';l150_withCorrection';l400_withCorrection';l100_withCorrection';l200_withCorrection';l50_withCorrection';l36_withCorrection';l25_withCorrection';l10_withCorrection'];


%pos = [60*ones(length(l400),1);90*ones(length(l100),1);4*60*ones(length(l50),1)];
pos_noCorrection = [120*ones(length(l400_noCorrection),1);120*ones(length(l400_noCorrection),1);120*ones(length(l400_noCorrection),1);120*ones(length(l400_noCorrection),1);180*ones(length(l200_noCorrection),1);180*ones(length(l100_noCorrection),1);4*60*ones(length(l50_noCorrection),1);4*60*ones(length(l36_noCorrection),1);4*60*ones(length(l25_noCorrection),1);;4*60*ones(length(l10_noCorrection),1);];

%pos_withCorrection = [];
%if length(l400_withCorrection)>0 
%  pos_withCorrection = [pos_withCorrection;auxData.numProcessors{2}(400)*ones(length(l400_withCorrection),1);];
%end
%if length(l200_withCorrection)>0
%  pos_withCorrection = [pos_withCorrection;auxData.numProcessors{2}(200)*ones(length(l200_withCorrection),1);];
%end
%if length(l100_withCorrection)>0
%  pos_withCorrection = [pos_withCorrection;auxData.numProcessors{2}(100)*ones(length(l100_withCorrection),1);];
%end
%if length(l50_withCorrection)>0
%  pos_withCorrection = [pos_withCorrection;auxData.numProcessors{2}(50)*ones(length(l50_withCorrection),1);];
%end
%if length(l25_withCorrection)>0
%  pos_withCorrection = [pos_withCorrection;auxData.numProcessors{2}(25)*ones(length(l25_withCorrection),1);];
%end

%keyboard

dup_noCorrection = [];
for i=1:length(l_noCorrection)-1
  for j=i+1:length(l_noCorrection)
    if length(l_noCorrection{i})==length(l_noCorrection{j})
      if strcmp(l_noCorrection{i},l_noCorrection{j})
        dup_noCorrection = [dup_noCorrection,j];
      end
    end
  end
end
l_noCorrection(dup_noCorrection) = []; 
pos_noCorrection(dup_noCorrection) = [];




name = ['tr1_2_',auxData.outFileDirName,'_noCorrection.txt']
fid = fopen([userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/',name],'w');
fprintf(fid,'#!/bin/bash \n');

k=1;
for i=1:length(l_noCorrection)
  w = [l_noCorrection{i},' ;sleep 20;'];
  fprintf(fid,'%s\n',w);
  if mod(k,5)==0
    fprintf(fid,'sleep %s;\n',[num2str(pos_noCorrection(i)),'m']);
  end
  k = k+1;
end
fclose(fid);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/',name])

%%%%%%%%
% with correction
% sum the number of CPU
name = ['tr1_2_',auxData.outFileDirName,'_withCorrection.txt']
fid = fopen([userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/',name],'w');
fprintf(fid,'#!/bin/bash \n');


k=1;
for i=1:length(l_withCorrection)
  w = [l_withCorrection{i},' ;sleep 20;'];
  fprintf(fid,'%s\n',w);
end
fclose(fid);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,'/',name])


