%%%%%%%%%
clear
readLength = 50;
basicSeqKey = '/homes/csfaculty/shental/CS/BAC/dat450000/key450000';
basicSeqNameDir = '/homes/csfaculty/shental/CS/BAC/dat450000/';
load('~/CS/BAC/dataForClose500','ind_bac_in_mix','correctWeight')

inifiniteNumberOfReadsFlag = 1;
[uniqueReads,uniqueReads_length,auxData.fracRelevantReadsForInfinity]=createReadsForInfiniteNumber(ind_bac_in_mix,correctWeight,readLength,basicSeqNameDir,basicSeqKey);
dataIn.fracRelevantReadsForInfinity = auxData.fracRelevantReadsForInfinity;

currData = [1:455055]';
currData(ind_bac_in_mix) = [];
currData = currData(randperm(length(currData)));

for i=155:length(ind_bac_in_mix)
  i
  currTmp= [ind_bac_in_mix(i);currData(1:999)];
  [normalizedBac values] = prepareGroupOf1000DistributedSequenceFiles2(readLength,currTmp,basicSeqNameDir,basicSeqKey);  
  [fracRelevantReads,sumRelevantReads(i)]=currReads(uniqueReads,uniqueReads_length,values,inifiniteNumberOfReadsFlag,dataIn);
  sumRelevantReads(i) = sum(fracRelevantReads);
  %fracRelevantReads = fracRelevantReads./10/sum(fracRelevantReads);
  [currResNoNorm{i},firstFlag]=runOneGroupOf1000post(1,normalizedBac,fracRelevantReads,0);
  A=full(normalizedBac'*normalizedBac);
  %D=rref_mod2(A);
  rnk(i) = rank(A);

  fracRelevantReads = fracRelevantReads./sum(fracRelevantReads);
  [currResNorm{i},firstFlag]=runOneGroupOf1000post(1,normalizedBac,fracRelevantReads,0);  
end


for i=1:50
  plot(currResNorm{i},currResNoNorm{i},'.')
  line([0,0.5],[0 0.5])
  pause
end


for i=1:50
  plot(currResNoNorm{i},'.')
  title(sumRelevantReads(i))
  pause
  
end

b = normalizedBac>0;
b = full(b);
yb = zeros(1000,1);
for j=1:1000
  yb(j) = sum(fracRelevantReads(find(b(:,j))));
end

z = b-b(:,1)*ones(1,1000);
zz = sum(abs(z),1);


plot(zz,ybx,'.')

% are these ok? - several are part of the correct 500 I think

[normalizedBacFull values] = prepareGroupOf1000DistributedSequenceFiles2(readLength,[ind_bac_in_mix;currData(1:999)],basicSeqNameDir,basicSeqKey);  
[fracRelevantReads,sumRelevantReads]=currReads(uniqueReads,uniqueReads_length,values,inifiniteNumberOfReadsFlag,dataIn);
b = normalizedBacFull>0;
b = full(b);
yb = zeros(size(normalizedBacFull,2),1);
for j=1:size(normalizedBacFull,2)
  yb(j) = sum(fracRelevantReads(find(b(:,j))));
end

zz = zeros(500,size(b,2));
for i=1:500
  i
  for j=1:size(b,2)
    zz(i,j) = sum(abs(b(:,i)-b(:,j)),1);
  end
end

minz = min(zz,[],1);

for i=1:50
  plot(currResNoNorm{i},minz([i,501:end]),'.')
  pause
end

% look for fn

cm = cell2mat(currResNoNorm);

for i=1:478
  if max(cm(:,i))>cm(1,i)
    i
    pause
  end
end

find(cm(1,:)<10^-2)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for Amnon:
% mix of close 500: 
clear
readLength = 50;
basicSeqKey = '/homes/csfaculty/shental/CS/BAC/dat450000/key450000';
basicSeqNameDir = '/homes/csfaculty/shental/CS/BAC/dat450000/';
load('~/CS/BAC/dataForClose500','ind_bac_in_mix','correctWeight')

inifiniteNumberOfReadsFlag = 1;
[uniqueReads,uniqueReads_length,auxData.fracRelevantReadsForInfinity]=createReadsForInfiniteNumber(ind_bac_in_mix,correctWeight,readLength,basicSeqNameDir,basicSeqKey);
dataIn.fracRelevantReadsForInfinity = auxData.fracRelevantReadsForInfinity;

currData = [1:455055]';
currData(ind_bac_in_mix) = [];
currData = currData(randperm(length(currData)));

% question which are supposed to be high - are they high since they are close to one of the correct bacteria?


[normalizedBac1 values1] = prepareGroupOf1000DistributedSequenceFiles2(readLength,[ind_bac_in_mix,;currData(1:1000)],basicSeqNameDir,basicSeqKey);  
inifiniteNumberOfReadsFlag = 1;
[fracRelevantReads1Infinite,sumRelevantReads]=currReads(uniqueReads,uniqueReads_length,values1,inifiniteNumberOfReadsFlag,dataIn);
[normalizedBac2 values2] = prepareGroupOf1000DistributedSequenceFiles2(readLength,[ind_bac_in_mix,;currData(1001:2000)],basicSeqNameDir,basicSeqKey);  
[fracRelevantReads2Infinite,sumRelevantReads]=currReads(uniqueReads,uniqueReads_length,values2,inifiniteNumberOfReadsFlag,dataIn);


%%%%%%%%%5
% finite number of reads
userDir = '~/';
auxData.reads_data_fileName = '~/CS/BAC/dataForSim/simSec500_Nbacmix_500_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads.mat';
[junk,red] = newReadsBasedSeed2(auxData.reads_data_fileName,userDir,basicSeqNameDir,basicSeqKey,500);
red = uint8(red);
[uniqueReads,uniqueReads_inds] = myUniqueUINT8(red);
clear red

uniqueReads_length = zeros(size(uniqueReads,1),1);
for i=1:size(uniqueReads_length,1)
  uniqueReads_length(i) = length(uniqueReads_inds{i});
end
clear uniqueReads_inds

dataIn = struct;
dataIn.normalizeByTotalNumberOfReadsFlag = 1;
[fracRelevantReads1Finite,sumRelevantReads]=currReads(uniqueReads,uniqueReads_length,values1,0,dataIn);
[fracRelevantReads2Finite,sumRelevantReads]=currReads(uniqueReads,uniqueReads_length,values2,0,dataIn);

save ~/CS/BAC/dataForAmnonInfiniteAndFinite normalizedBac1 normalizedBac2 fracRelevantReads1Finite fracRelevantReads2Finite fracRelevantReads1Infinite fracRelevantReads2Infinite correctWeight 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finite number of reads
% mix of close 500: 
clear
readLength = 50;
basicSeqKey = '/homes/csfaculty/shental/CS/BAC/dat450000/key450000';
basicSeqNameDir = '/homes/csfaculty/shental/CS/BAC/dat450000/';
load('~/CS/BAC/dataForClose500','ind_bac_in_mix','correctWeight')

currData = [1:455055]';
currData(ind_bac_in_mix) = [];
currData = currData(randperm(length(currData)));

[normalizedBac3 values1] = prepareGroupOf1000DistributedSequenceFiles2(readLength,[ind_bac_in_mix,;currData(1:1000)],basicSeqNameDir,basicSeqKey);  

userDir = '~/';


auxData.reads_data_fileName = '~/CS/BAC/dataForSim/simSec500_Nbacmix_500_Nread_100000_Readlen_50_Npower_05_bacdistflag_1_NoReads.mat';
[junk,red] = newReadsBasedSeed2(auxData.reads_data_fileName,userDir,basicSeqNameDir,basicSeqKey,500);
red = uint8(red);
[uniqueReads,uniqueReads_inds] = myUniqueUINT8(red);
clear red

uniqueReads_length = zeros(size(uniqueReads,1),1);
for i=1:size(uniqueReads_length,1)
  uniqueReads_length(i) = length(uniqueReads_inds{i});
end
clear uniqueReads_inds

dataIn = struct;
dataIn.normalizeByTotalNumberOfReadsFlag = 1;
[fracRelevantReadsFinite100000,sumRelevantReads]=currReads(uniqueReads,uniqueReads_length,values1,0,dataIn);

auxData.reads_data_fileName = '~/CS/BAC/dataForSim/simSec500_Nbacmix_500_Nread_20000_Readlen_50_Npower_05_bacdistflag_1_NoReads.mat';
[junk,red] = newReadsBasedSeed2(auxData.reads_data_fileName,userDir,basicSeqNameDir,basicSeqKey,500);
red = uint8(red);
[uniqueReads,uniqueReads_inds] = myUniqueUINT8(red);
clear red

uniqueReads_length = zeros(size(uniqueReads,1),1);
for i=1:size(uniqueReads_length,1)
  uniqueReads_length(i) = length(uniqueReads_inds{i});
end
clear uniqueReads_inds

dataIn = struct;
dataIn.normalizeByTotalNumberOfReadsFlag = 1;
[fracRelevantReadsFinite20000,sumRelevantReads]=currReads(uniqueReads,uniqueReads_length,values1,0,dataIn);


save ~/CS/BAC/dataForAmnonFinite normalizedBac3 fracRelevantReadsFinite100000 fracRelevantReadsFinite20000


% for Amnon:
% mix of distant 500: 
clear
readLength = 50;
basicSeqKey = '/homes/csfaculty/shental/CS/BAC/dat450000/key450000';
basicSeqNameDir = '/homes/csfaculty/shental/CS/BAC/dat450000/';
load('~/CS/BAC/dataForClose500','ind_bac_in_mix','correctWeight')

rand('state',1047)
ind_bac_in_mix = randperm(455055);
ind_bac_in_mix = ind_bac_in_mix(1:500)';

c = zeros(size(correctWeight));
c(ind_bac_in_mix) = correctWeight(find(correctWeight));
correctWeight = c;

freq = [1:500;correctWeight(ind_bac_in_mix)'];
save ~/CS/BAC/freq freq

% save ~/CS/mFilesBAC/ii correctWeight ind_bac_in_mix 
% load ~/CS/BAC/dataForSim/simSec500_Nbacmix_500_Nread_2000000_Readlen_50_Npower_05_bacdistflag_0_NoReads;clear w;load ~/CS/mFilesBAC/ii;save ~/CS/BAC/dataForSim/simSec500_Nbacmix_500_Nread_2000000_Readlen_50_Npower_05_bacdistflag_0_NoReads

inifiniteNumberOfReadsFlag = 1;
[uniqueReads,uniqueReads_length,auxData.fracRelevantReadsForInfinity]=createReadsForInfiniteNumber(ind_bac_in_mix,correctWeight,readLength,basicSeqNameDir,basicSeqKey);
dataIn.fracRelevantReadsForInfinity = auxData.fracRelevantReadsForInfinity;

currData = [1:455055]';
currData(ind_bac_in_mix) = [];
currData = currData(randperm(length(currData)));

% question which are supposed to be high - are they high since they are close to one of the correct bacteria?


[normalizedBac1 values1] = prepareGroupOf1000DistributedSequenceFiles2(readLength,[ind_bac_in_mix;currData(1:1000)],basicSeqNameDir,basicSeqKey);  
inifiniteNumberOfReadsFlag = 1;
[fracRelevantReads1Infinite,sumRelevantReads]=currReads(uniqueReads,uniqueReads_length,values1,inifiniteNumberOfReadsFlag,dataIn);
[normalizedBac2 values2] = prepareGroupOf1000DistributedSequenceFiles2(readLength,[ind_bac_in_mix,;currData(1001:2000)],basicSeqNameDir,basicSeqKey);  
[fracRelevantReads2Infinite,sumRelevantReads]=currReads(uniqueReads,uniqueReads_length,values2,inifiniteNumberOfReadsFlag,dataIn);

%%%%%%%%%5
% finite number of reads
userDir = '~/';
auxData.reads_data_fileName = '~/CS/BAC/dataForSim/simSec500_Nbacmix_500_Nread_2000000_Readlen_50_Npower_05_bacdistflag_0_NoReads.mat';

[junk,red] = newReadsBasedSeed2(auxData.reads_data_fileName,userDir,basicSeqNameDir,basicSeqKey,500);
red = uint8(red);
[uniqueReads,uniqueReads_inds] = myUniqueUINT8(red);
clear red

uniqueReads_length = zeros(size(uniqueReads,1),1);
for i=1:size(uniqueReads_length,1)
  uniqueReads_length(i) = length(uniqueReads_inds{i});
end
clear uniqueReads_inds

dataIn = struct;
dataIn.normalizeByTotalNumberOfReadsFlag = 1;
[fracRelevantReads1Finite,sumRelevantReads]=currReads(uniqueReads,uniqueReads_length,values1,0,dataIn);
[fracRelevantReads2Finite,sumRelevantReads]=currReads(uniqueReads,uniqueReads_length,values2,0,dataIn);

save ~/CS/BAC/dataForAmnonInfiniteAndFiniteDistant500 normalizedBac1 normalizedBac2 fracRelevantReads1Finite fracRelevantReads2Finite fracRelevantReads1Infinite fracRelevantReads2Infinite correctWeight 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% try random seq + follow

