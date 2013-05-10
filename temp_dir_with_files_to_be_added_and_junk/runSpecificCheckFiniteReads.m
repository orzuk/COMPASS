function runSpecificCheckFiniteReads(outFileDirName,numIter,readLength,pwr,nb,nr,thresholdForCollectingBAC,groupSize)

if ~exist('groupSize')
  groupSize = 1000;
end
groupSize
if ~exist('thresholdForCollectingBAC')
  thresholdForCollectingBAC = 10^-3;
end
thresholdForCollectingBAC

userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
addpath([userDir,'/CS/mFilesBAC'])
unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/',outFileDirName])
prefix = outFileDirName;
brFlag = 1;
crGeneralFourth(prefix,outFileDirName,10,[0],nb,pwr,numIter,nr,[0 1],[0 1],readLength,'hour','week',150000,50000,400,brFlag,1,thresholdForCollectingBAC,groupSize);
crGeneralFourth(prefix,outFileDirName,10,[0],nb,pwr,numIter,nr,[0 1],[0 1],readLength,'hour','week',150000,50000,400,brFlag,0,thresholdForCollectingBAC,groupSize);

[junk,list50] = unix(['grep -l "readlen_50" ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/list*']);
[junk,list100] = unix(['grep -l "readlen_100" ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/list*']);
[junk,list400] = unix(['grep -l "readlen_400" ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/list*']);

l50 = extract_listNames(list50);
l100 = extract_listNames(list100);
l400 = extract_listNames(list400);

l = [l400';l100';l50';];
%pos = [60*ones(length(l400),1);90*ones(length(l100),1);4*60*ones(length(l50),1)];
pos = [120*ones(length(l400),1);180*ones(length(l100),1);4*60*ones(length(l50),1)];

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
fid = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name],'w');
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
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name])
