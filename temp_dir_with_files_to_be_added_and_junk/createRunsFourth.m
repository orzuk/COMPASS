%%%%%%%%%%%%%%%%%%
clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
addpath([userDir,'/CS/mFilesBAC'])
outFileDirName = 'ffFourth';
unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/',outFileDirName])
prefix = outFileDirName;
brFlag = 1;

if 1==2
  unix(['rm -fr ',userDir,'/CS/tmp/',outFileDirName,';',...
        ';rm -fr ',userDir,'/CS/BAC/dataForSim/',outFileDirName,';',...
        ';rm -fr ',userDir,'/CS/BAC/tmpRuns/',outFileDirName,';',...
        ';rm -fr ',userDir,'/CS/BAC/',outFileDirName,';',...
        ';rmdir ',userDir,'/CS/tmp/',outFileDirName,';',...
        ';rmdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,';',...
        ';rmdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName,';',...
        ';rmdir ',userDir,'/CS/BAC/',outFileDirName,';'])
end


crGeneralFourth(prefix,outFileDirName,10,[0],[10 100 1000],[0 0.5],20,[inf],0,0,[50 100 400],'hour','week',150000,20000,400,brFlag,1);
crGeneralFourth(prefix,outFileDirName,10,[0],[10 100 1000],[0 0.5],20,[inf],0,0,[50 100 400],'hour','week',150000,20000,400,brFlag,0);

[junk,list50] = unix(['grep -l "readlen_50" ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/list*']);
[junk,list100] = unix(['grep -l "readlen_100" ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/list*']);
[junk,list400] = unix(['grep -l "readlen_400" ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/list*']);
[junk,list700] = unix(['grep -l "readlen_700" ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/list*']);

% list them to find intersect
l50 = extract_listNames(list50);
l100 = extract_listNames(list100);
l400 = extract_listNames(list400);
l700 = extract_listNames(list700);

l = [l50';l100';l400';l700'];
pos = [120*ones(length(l50),1);60*ones(length(l100),1);40*ones(length(l400),1);30*ones(length(l700),1)];
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

outFileDirName = 'ffFourth';
userDir = getuserdir;
%userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
basicDir = [userDir,'/CS/BAC/',outFileDirName,'/'];

clear res
nb = 10;
for readLength=[50 100 400]
  for np=[0 0.5]
    res(readLength,np*10+1,:) = extractRes(basicDir,outFileDirName,0,nb,np,Inf,0,0,readLength,410849,1);
  end
end

nb = 10
%%%%%%%%%%%%%%%%%%%
clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
addpath([userDir,'/CS/mFilesBAC'])
outFileDirName = 'ffFourth2';
unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/',outFileDirName])
prefix = outFileDirName;
brFlag = 1;

if 1==2
  unix(['rm -fr ',userDir,'/CS/tmp/',outFileDirName,';',...
        ';rm -fr ',userDir,'/CS/BAC/dataForSim/',outFileDirName,';',...
        ';rm -fr ',userDir,'/CS/BAC/tmpRuns/',outFileDirName,';',...
        ';rm -fr ',userDir,'/CS/BAC/',outFileDirName,';',...
        ';rmdir ',userDir,'/CS/tmp/',outFileDirName,';',...
        ';rmdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,';',...
        ';rmdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName,';',...
        ';rmdir ',userDir,'/CS/BAC/',outFileDirName,';'])
end

crGeneralFourth(prefix,outFileDirName,10,[0],[10 100 1000],[0 0.5],80,[inf],0,0,[50 100 400],'hour','week',150000,20000,400,brFlag,1);
crGeneralFourth(prefix,outFileDirName,10,[0],[10 100 1000],[0 0.5],80,[inf],0,0,[50 100 400],'hour','week',150000,20000,400,brFlag,0);

[junk,list50] = unix(['grep -l "readlen_50" ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/list*']);
[junk,list100] = unix(['grep -l "readlen_100" ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/list*']);
[junk,list400] = unix(['grep -l "readlen_400" ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/list*']);
[junk,list700] = unix(['grep -l "readlen_700" ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/list*']);

% list them to find intersect
l50 = extract_listNames(list50);
l100 = extract_listNames(list100);
l400 = extract_listNames(list400);
l700 = extract_listNames(list700);

l = [l50';l100';l400';l700'];
pos = [240*ones(length(l50),1);120*ones(length(l100),1);80*ones(length(l400),1);60*ones(length(l700),1)];
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


%%%%%%%%%%%%%%%
outFileDirName = 'ffFourth2';
userDir = getuserdir;
%userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
basicDir = [userDir,'/CS/BAC/',outFileDirName,'/'];

clear res
nb = 10;
for readLength=[50 100 400]
  for np=[0 0.5]
    res(readLength,np*10+1,:) = extractRes(basicDir,outFileDirName,0,nb,np,Inf,0,0,readLength,410849,1);
  end
end

%%%%%%%%%%%%%%%555
% run50 with higher number or repeats
clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
addpath([userDir,'/CS/mFilesBAC'])
outFileDirName = 'ffJust50';
unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/',outFileDirName])
prefix = outFileDirName;
brFlag = 1;

crGeneralFourth(prefix,outFileDirName,10,[0],[100 1000],[0 0.5],100,[inf],0,0,[50],'hour','week',150000,50000,400,brFlag,1);
crGeneralFourth(prefix,outFileDirName,10,[0],[100 1000],[0 0.5],100,[inf],0,0,[50],'hour','week',150000,50000,400,brFlag,0);

[junk,list50] = unix(['grep -l "readlen_50" ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/list*']);

% list them to find intersect
l50 = extract_listNames(list50);


name = 'tr1_2.txt'
fid = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name],'w');
fprintf(fid,'#!/bin/bash \n');

k=1;
for i=1:length(l50)
  w = [l50{i},' ;sleep 20;'];
  fprintf(fid,'%s\n',w);
  if mod(k,5)==0
    fprintf(fid,'sleep 240m\n');
  end
  k = k+1;
end
fclose(fid);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name])

userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
basicDir = [userDir,'/CS/BAC/',outFileDirName,'/'];

clear res
nb = 100;
for readLength=[50]
  for np=[0 0.5]
    res(readLength,np*10+1,:) = extractRes(basicDir,outFileDirName,0,nb,np,Inf,0,0,readLength,410849,1);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%

% run50 with higher number or repeats
4:30 in the morning

clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
addpath([userDir,'/CS/mFilesBAC'])
outFileDirName = 'newKmer20';
unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/',outFileDirName])
prefix = outFileDirName;
brFlag = 1;

if 1==2
  unix(['rm -fr ',userDir,'/CS/tmp/',outFileDirName,';',...
        ';rm -fr ',userDir,'/CS/BAC/dataForSim/',outFileDirName,';',...
        ';rm -fr ',userDir,'/CS/BAC/tmpRuns/',outFileDirName,';',...
        ';rm -fr ',userDir,'/CS/BAC/',outFileDirName,';',...
        ';rmdir ',userDir,'/CS/tmp/',outFileDirName,';',...
        ';rmdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,';',...
        ';rmdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName,';',...
        ';rmdir ',userDir,'/CS/BAC/',outFileDirName,';'])
end


crGeneralFourth(prefix,outFileDirName,10,[0],[100],[0 0.5],20,[inf],0,0,[50 100 400],'hour','week',150000,50000,400,brFlag,1);
crGeneralFourth(prefix,outFileDirName,10,[0],[100],[0 0.5],20,[inf],0,0,[50 100 400],'hour','week',150000,50000,400,brFlag,0);

[junk,list50] = unix(['grep -l "readlen_50" ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/list*']);
[junk,list100] = unix(['grep -l "readlen_100" ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/list*']);
[junk,list400] = unix(['grep -l "readlen_400" ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/list*']);

l50 = extract_listNames(list50);
l100 = extract_listNames(list100);
l400 = extract_listNames(list400);

l = [l400';l100';l50';];
pos = [60*ones(length(l400),1);90*ones(length(l100),1);4*60*ones(length(l50),1)];
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




userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
basicDir = [userDir,'/CS/BAC/',outFileDirName,'/'];

clear res
nb = 100;
for readLength=[50 100 400]
  for np=[0 0.5]
    res(readLength,np*10+1,:) = extractRes(basicDir,outFileDirName,0,nb,np,Inf,0,0,readLength,410849,1);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%55
clear

runSpecific('test05LongReads',100,[100 400],[0.5],100)

%1000 with 10-4
runSpecific('test05LongReads_numBacteria_1000',1,[400],[0],1000,10^-4,5000)
% got stuck

%%%%%%%%%%%%%
% test finite reads with/without correction 
runSpecificCheckFiniteReads('testFiniteReads',1,[50],[0 0.5],[10 100 1000],[10^4 10^5 10^6])

runSpecificCheckFiniteReads('testFiniteReads2',10,[50],[0 0.5],[10 100 1000],[10^4 10^5 10^6])

outFileDirName = 'testFiniteReads2'
userDir = getuserdir;
%userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
basicDir = [userDir,'/CS/BAC/',outFileDirName,'/'];


clear res
nb = [10];
for readLength=[50]
  for np=[0]%[0 0.5]
    for nr=[10^5]%[10^4 10^5 10^6]
      for addNoiseFlagInput=[0 1]
        for correctMeasurementForReadErrorsFlagInput=[0,1]
          if ~(addNoiseFlagInput==0 & correctMeasurementForReadErrorsFlagInput==1)
            res(readLength,np*10+1,:) = extractRes(basicDir,outFileDirName,0,nb,np,nr,addNoiseFlagInput,correctMeasurementForReadErrorsFlagInput,readLength,410849,1);
          end
        end
      end
    end
  end
end


% check close/far


%%%%%%%%%%%%%%%%%%%
%%%%%%%%%
% close bacteria reads - infinite reads
auxData = struct;
auxData.outFileDirName = 'close';
auxData.bdf = 1;
auxData.nb = [10 100 1000];
auxData.pwr = [0 0.5];
auxData.readLength = [50 100 400];
auxData.numIter = 2;
auxData.nr = [Inf];
auxData.addNoiseFlagInput = 0;
auxData.correctMeasurementForReadErrorsFlagInput = 0;

runGeneral(auxData);


%%%%%%%%%%
% distant bacteria - infinite number of reads, short read length
auxData = struct;
auxData.outFileDirName = 'shortReadLengthInfReads';
auxData.bdf = 0;
auxData.nb = [10 100 1000];
auxData.pwr = [0 0.5];
auxData.readLength = [36];
auxData.numIter = 20;
auxData.nr = [Inf];
auxData.addNoiseFlagInput = 0;
auxData.correctMeasurementForReadErrorsFlagInput = 0;

runGeneral(auxData);


clear res
nb = [100];
for readLength=[36]
  for np=[0 0.5]
    for nr=[Inf]
      for addNoiseFlagInput=[0]
        for correctMeasurementForReadErrorsFlagInput=[0]
          if ~(addNoiseFlagInput==0 & correctMeasurementForReadErrorsFlagInput==1)
            res(readLength,np*10+1,:) = extractRes(basicDir,outFileDirName,0,nb,np,nr,addNoiseFlagInput,correctMeasurementForReadErrorsFlagInput,readLength,410849,1);
          end
        end
      end
    end
  end
end


% compare results with/without noise


% distant bacteria - infinite number of reads, short read length
auxData = struct;
auxData.outFileDirName = 'shortReadLengthInfReads25';
auxData.bdf = 0;
auxData.nb = [10 100 1000];
auxData.pwr = [0 0.5];
auxData.readLength = [25];
auxData.numIter = 10;
auxData.nr = [Inf];
auxData.addNoiseFlagInput = 0;
auxData.correctMeasurementForReadErrorsFlagInput = 0;

runGeneral(auxData);

dispExtractRes(auxData,'~/',10)
dispExtractRes(auxData,'~/',100)
dispExtractRes(auxData,'~/',1000)

%%%%%%%%%%%%%%%%%%%%%%%%
clear 
auxData = struct;
auxData.outFileDirName = 'shortReadLengthInfReads10';
auxData.bdf = 0;
auxData.nb = [10 100 1000];
auxData.pwr = [0 0.5];
auxData.readLength = [10];
auxData.numIter = 10;
auxData.nr = [Inf];
auxData.addNoiseFlagInput = 0;
auxData.correctMeasurementForReadErrorsFlagInput = 0;

runGeneral(auxData);



%%%%%%%%%%%%%%%%%%55
% distant bacteria - find  number of reads, long read length
clear
auxData = struct;
auxData.outFileDirName = 'finiteReadsLong';
auxData.bdf = 0;
auxData.nb = [100 1000];
auxData.pwr = [0 0.5];
auxData.readLength = [100 400];
auxData.numIter = 10;
auxData.nr = [10^4 10^5 10^6];
auxData.addNoiseFlagInput = [0 1];
auxData.correctMeasurementForReadErrorsFlagInput = 0;

runGeneral(auxData);


%%%%%%%%%%%%%%%%%%55
% distant bacteria - find  number of reads, read length 100
clear
auxData = struct;
auxData.outFileDirName = 'finiteReadsLong400Second';
auxData.bdf = 0;
auxData.nb = [100 1000];
auxData.pwr = [0 0.5];
auxData.readLength = [400];
auxData.numIter = 1;
auxData.nr = [10^4 10^5 10^6];
auxData.addNoiseFlagInput = [0 1];
auxData.correctMeasurementForReadErrorsFlagInput = 0;

runGeneral(auxData);

% distant bacteria - find  number of reads, read length 100
clear
auxData = struct;
auxData.outFileDirName = 'finiteReadsLong400Third';
auxData.bdf = 0;
auxData.nb = [100 1000];
auxData.pwr = [0 0.5];
auxData.readLength = [400];
auxData.numIter = 20;
auxData.nr = [10^4 10^5];
auxData.addNoiseFlagInput = [0 1];
auxData.correctMeasurementForReadErrorsFlagInput = 0;

runGeneral(auxData);


% distant bacteria - find  number of reads, read length 100
clear
auxData = struct;
auxData.outFileDirName = 'finiteReadsLong400for10e6';
auxData.bdf = 0;
auxData.nb = [100 1000];
auxData.pwr = [0 0.5];
auxData.readLength = [400];
auxData.numIter = 20;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [0 1];
auxData.correctMeasurementForReadErrorsFlagInput = 0;

runGeneral(auxData);

%%
dd
userDir = getuserdir;
auxData.basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];
auxData.basicSeqKey = [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];

dr = dir(['~/CS/BAC/',auxData.outFileDirName,'/',auxData.outFileDirName,'_*']);
matlabpool open 6
parfor (i=1:length(dr))
 CreateMothurDistNoam(['~/CS/BAC/',auxData.outFileDirName,'/'],dr(i).name,['mothurRes_',auxData.outFileDirName,'/mothur_',dr(i).name],0,auxData);
end



%%%%%%%%%%%%%%5
% finite reads 400 10^6

clear
auxData = struct;
auxData.outFileDirName = 'finiteReadsLong400for10e6';
auxData.bdf = 0;
auxData.nb = [100 1000];
auxData.pwr = [0 0.5];
auxData.readLength = [400];
auxData.numIter = 20;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [0 1];
auxData.correctMeasurementForReadErrorsFlagInput = 0;

runGeneral(auxData);

% add reads to 400 and 10^6
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
doLocalCreateReadsForSpecificCases(auxData,userDir)
% ./tr1


basicDir = [userDir,'/CS/BAC/',auxData.outFileDirName,'/'];

extractRes(basicDir,auxData.outFileDirName,0,100,0.5,1000000,0,0,400,410849,1);


%%%%%%%%%%%%%%5
% finite reads trySetting

clear
auxData = struct;
auxData.outFileDirName = 'trySetting';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [0 0.5];
auxData.readLength = [100];
auxData.numIter = 20;
auxData.nr = [10^5 5*10^5 10^6];
auxData.addNoiseFlagInput = [0 1];
auxData.correctMeasurementForReadErrorsFlagInput = 0;

runGeneral(auxData);

% add reads to 400 and 10^6
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
doLocalCreateReadsForSpecificCases(auxData,userDir)
% ./tr1

%
userDir = getuserdir;
auxData.basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];
auxData.basicSeqKey = [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];

dr = dir(['~/CS/BAC/',auxData.outFileDirName,'/',auxData.outFileDirName,'_*']);
matlabpool open 4
parfor (i=1:length(dr))
 CreateMothurDistNoam(['~/CS/BAC/',auxData.outFileDirName,'/'],dr(i).name,['mothurRes_',auxData.outFileDirName,'/mothur_',dr(i).name],0,auxData);
end



%%%%%%%%%%%555
clear
auxData = struct;
auxData.outFileDirName = 'infinite1000';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [0 0.5];
auxData.readLength = [50 100];
auxData.numIter = 20;
auxData.nr = [inf];
auxData.addNoiseFlagInput = [0];
auxData.correctMeasurementForReadErrorsFlagInput = 0;

runGeneral(auxData);

% add reads to 400 and 10^6
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
doLocalCreateReadsForSpecificCases(auxData,userDir)
% ./tr1

userDir = getuserdir;
auxData.basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];
auxData.basicSeqKey = [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
dr = dir(['~/CS/BAC/',auxData.outFileDirName,'/',auxData.outFileDirName,'_*']);

matlabpool open 6

parfor (i=1:length(dr))
 CreateMothurDistNoam(['~/CS/BAC/',auxData.outFileDirName,'/'],dr(i).name,['mothurRes_',auxData.outFileDirName,'/mothur_',dr(i).name],0,auxData);
end
matlabpool close


%%%%%%%%%%%%%%%%%%%%%%%%%55
clear

auxData = struct;
auxData.brFlag = 0;
auxData.outFileDirName = 'testFiniteReads2';

userDir = getuserdir;
%dr = dir([userDir,'/CS/BAC/',auxData.outFileDirName,'/*Nbac_in_mixture_100_*Nreads_10000*']);
dr = dir([userDir,'/CS/BAC/',auxData.outFileDirName,'/*Nbac_in_mixture_100_*Nreads_10000010000*']);



auxData.numProcessors = 1;
auxData.dr = dr;

%distribute_Mah(auxData)

clear
userDir = getuserdir;

auxData.outFileDirName = 'testFiniteReads2';
auxData.basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];
auxData.basicSeqKey = [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];

dr = dir(['~/CS/BAC/',auxData.outFileDirName,'/',auxData.outFileDirName,'_*']);
for i=1:length(dr)
  %currName = dr(i).name;
  %a = findstr(currName,'.mat');
  
  CreateMothurDistNoam(['~/CS/BAC/',auxData.outFileDirName,'/'],dr(i).name,['mothurRes_',auxData.outFileDirName,'/mothur_',dr(i).name],0,auxData);
end


matlabpool open 4
parfor (i=82:length(dr))
 CreateMothurDistNoam(['~/CS/BAC/',auxData.outFileDirName,'/'],dr(i).name,['mothurRes_',auxData.outFileDirName,'/mothur_',dr(i).name],0,auxData);
end

% testing


ResultsDistance('~/CS/BAC/testFiniteReads2/','testFiniteReads2_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_9_Nreads_10000_Noise_1_Correction_0',auxData);

currNumProcessors = 4;
n = length(auxData.dr);
part = 1:currNumProcessors:n;
part(end) = n+1;
tmpInd = cell(length(part)-1,1);
for i=1:length(part)-1
  tmpInd{i} = part(i):part(i+1)-1;
end


matlabpool open local 4
parfor i=1:length(tmpInd)
%for i=1:length(tmpInd)
  for j=tmpInd{i}
    ResultsDistance([userDir,'/CS/BAC/',auxData.outFileDirName,'/'],auxData.dr(j).name,auxData);
  end  
end



%%%%%%%%%%%
clear
auxData = struct;
auxData.outFileDirName = 'test1';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [0];
auxData.readLength = [100];
auxData.numIter = 1;
auxData.nr = [inf];
auxData.addNoiseFlagInput = [0];
auxData.correctMeasurementForReadErrorsFlagInput = 0;

runGeneral(auxData);

% add reads to 400 and 10^6
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
doLocalCreateReadsForSpecificCases(auxData,userDir)
% ./tr1

%%%%%%%%%%%%%%%%%%%%
% test primers750
clear
auxData = struct;
auxData.outFileDirName = 'p750';
auxData.bdf = 0;
auxData.nb = [100 1000];
auxData.pwr = [0 0.5];
auxData.readLength = [50 100 400];
auxData.numIter = 1;
auxData.nr = [10^4 10^5 5*10^5 10^6 inf];
auxData.addNoiseFlagInput = [0 1];
auxData.correctMeasurementForReadErrorsFlagInput = 0;
auxData.dataSetPrimers750Flag = 1;
runGeneral(auxData);

% add reads to 400 and 10^6
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
doLocalCreateReadsForSpecificCases(auxData,userDir)
% ./tr1

% test primers750
clear
auxData = struct;
auxData.outFileDirName = 'p750Sec';
auxData.bdf = 0;
auxData.nb = [100 1000];
auxData.pwr = [0 0.5];
auxData.readLength = [50 100 400];
auxData.numIter = 10;
auxData.nr = [10^4 10^5 5*10^5 10^6 inf];
auxData.addNoiseFlagInput = [0 1];
auxData.correctMeasurementForReadErrorsFlagInput = 0;
auxData.dataSetPrimers750Flag = 1;
runGeneral(auxData);





% add reads to 400 and 10^6
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
doLocalCreateReadsForSpecificCases(auxData,userDir)
% ./tr1

 bjobs -u orzuk | grep hour | grep PEND|wc|cut -f1
 
 
 
distributeBAC_generalOrFourth('/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/p750Sec/structFor_p750Sec_bac_dist_flag_0_Nbac_in_mixture_1000_npower_5_readlen_50_numIter_3_Nreads_500000_Noise_1_Correction_0')


%%%%%%%%%%%%%%%%%%%%%%%%
% test new correction - method 1
clear
auxData = struct;
auxData.outFileDirName = 'testCorrection_method1';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [0];
auxData.readLength = [50];
auxData.numIter = 1;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = 1;
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors = 20;
runGeneral(auxData);
% end method 1






% 400
clear
auxData = struct;
auxData.outFileDirName = 'testCorrection_Ver1_readLen400';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [0];
auxData.readLength = [400];
auxData.numIter = 1;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = 1;
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors = 10;
runGeneral(auxData);

userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
doLocalCreateReadsForSpecificCases(auxData,userDir)

clear
auxData = struct;
auxData.outFileDirName = 'testCorrection_Ver2_readLen400';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [0];
auxData.readLength = [400];
auxData.numIter = 1;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = 1;
auxData.CorrectReadsNew = 2; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors = 10;
runGeneral(auxData);

userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
doLocalCreateReadsForSpecificCases(auxData,userDir)



% ran 50 with/without correction
clear
auxData = struct;
auxData.outFileDirName = 'testCorrection_method50';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [0];
auxData.readLength = [50];
auxData.numIter = 1;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [0 1];
auxData.correctMeasurementForReadErrorsFlagInput = [0 1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors = 20;
runGeneral(auxData);
%%%%%%%%%%%%%%% done


% ran 100 with/without correction
clear
auxData = struct;
auxData.outFileDirName = 'testCorrection_method100';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [0];
auxData.readLength = [100];
auxData.numIter = 1;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [0 1];
auxData.correctMeasurementForReadErrorsFlagInput = [0 1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors = [10 35];
runGeneral(auxData);

% ran 400 with/without correction
clear
auxData = struct;
auxData.outFileDirName = 'testCorrection_method400';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [0];
auxData.readLength = [400];
auxData.numIter = 1;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [0 1];
auxData.correctMeasurementForReadErrorsFlagInput = [0 1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors = [10 35]; % the second number refers to the correction case!
runGeneral(auxData);

userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
doLocalCreateReadsForSpecificCases(auxData,userDir)

