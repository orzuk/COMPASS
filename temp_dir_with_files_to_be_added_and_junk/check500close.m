clear

w = [];
%%%%%
% run with original order
auxData = struct;
  
auxData.tmpFileName = ['testSimSec500OriginalOrder'];
auxData.keepOriginalOrderFlag = 1;
auxData.saveName = ['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/',auxData.tmpFileName];

auxData.createReadsFlag = 1;
auxData.randomizeSeedAfterCreateReadsFlag = 0;
auxData.brFlag = 1;
auxData.queueName = 'hour';
auxData.reads_data_fileName = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/simSec500_Nbacmix_500_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads.mat';
auxData.readLength = 50;
auxData.numProcessors = 100;
auxData.firstFlag = 0;
auxData.groupSize = 1000;
auxData.smallestSetOfCollected = 1000;
auxData.thresholdForCollectingBAC = 1e-3;
auxData.numBACtoConsider = 455055;
auxDataFileName = ['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/structFor_',auxData.tmpFileName];
save(auxDataFileName,'auxData');
w = [w,';',sprintf(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/cvx/run_distributeBAC_general.sh  /broad/software/nonfree/Linux/' ...
                'redhat_5_x86_64/pkgs/matlab_2010b %s;'],auxDataFileName)];


%%%%%%%%%%%%%%%
% run with lower threshold

auxData = struct;
  
%%%
auxData.tmpFileName = ['testSimSec500LowerThreshold'];
auxData.thresholdForCollectingBAC = 1e-4;
auxData.randomizeSeedAfterCreateReadsFlag = 0;
auxData.saveName = ['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/',auxData.tmpFileName];

auxData.createReadsFlag = 1;
auxData.brFlag = 1;
auxData.queueName = 'hour';
auxData.reads_data_fileName = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/simSec500_Nbacmix_500_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads.mat';
auxData.readLength = 50;
auxData.numProcessors = 100;
auxData.firstFlag = 0;
auxData.groupSize = 1000;
auxData.smallestSetOfCollected = 1000;
auxData.numBACtoConsider = 455055;
auxDataFileName = ['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/structFor_',auxData.tmpFileName];
save(auxDataFileName,'auxData');
w = [w,';',sprintf(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/cvx/run_distributeBAC_general.sh  /broad/software/nonfree/Linux/' ...
                'redhat_5_x86_64/pkgs/matlab_2010b %s;'],auxDataFileName)];

%%%%%%%%%%%%%%%%%%%55
% 

% run 10 times with different seeds


for i=1:10
  auxData = struct;
  
  %%%
  auxData.tmpFileName = ['testSimSec500_round',num2str(i)];
  auxData.randomizeSeedAfterCreateReadsFlag = 1;
  auxData.saveName = ['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/',auxData.tmpFileName];

  %%
  
  auxData.createReadsFlag = 1;
  auxData.brFlag = 1;
  auxData.queueName = 'hour';
  auxData.reads_data_fileName = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/simSec500_Nbacmix_500_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads.mat';
  auxData.readLength = 50;
  auxData.numProcessors = 100;
  auxData.firstFlag = 0;
  auxData.groupSize = 1000;
  auxData.smallestSetOfCollected = 1000;
  auxData.numBACtoConsider = 455055;
  auxData.thresholdForCollectingBAC = 1e-3;
  auxDataFileName = ['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/structFor_',auxData.tmpFileName];
  save(auxDataFileName,'auxData');
  
  w = [w,sprintf(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/cvx/run_distributeBAC_general.sh  /broad/software/nonfree/Linux/' ...
                  'redhat_5_x86_64/pkgs/matlab_2010b %s;'],auxDataFileName)];


end


%%%%%%%% run with higher read length
ww = [];
for readLength=100:-10:60
  load('/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/simSec500_Nbacmix_500_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads.mat');
  readlen = readLength;
  fileName = ['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/simSec500_Nbacmix_500_Nread_2000000_Readlen_',num2str(readlen),'_Npower_05_bacdistflag_1_NoReads.mat']; 
  save(fileName,'Nbac_in_mixture','correctWeight','outseed','Nreads','ind_bac_in_mix','readlen','bac_dist_flag','npower','w');
                 
                      
                
  auxData = struct;
  
  %%%
  auxData.tmpFileName = ['testSimSec500_readLength',num2str(readLength)];
  auxData.randomizeSeedAfterCreateReadsFlag = 1;
  auxData.saveName = ['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/',auxData.tmpFileName];

  %%
  
  auxData.createReadsFlag = 1;
  auxData.brFlag = 1;
  auxData.queueName = 'hour';
  auxData.reads_data_fileName = fileName;
  auxData.readLength = readLength;
  auxData.numProcessors = 100;
  auxData.firstFlag = 0;
  auxData.groupSize = 1000;
  auxData.smallestSetOfCollected = 1000;
  auxData.numBACtoConsider = 455055;
  auxData.thresholdForCollectingBAC = 1e-3;
  auxDataFileName = ['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/structFor_',auxData.tmpFileName];
  save(auxDataFileName,'auxData');
  
  ww = [ww,sprintf(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/cvx/run_distributeBAC_general.sh  /broad/software/nonfree/Linux/' ...
                  'redhat_5_x86_64/pkgs/matlab_2010b %s;'],auxDataFileName)];


end





%%%%%%%%%%


plotRes('~/CS/BAC/dataForSim/simSec500_Nbacmix_500_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads.mat','~/CS/BAC/testSimSec500LowerThreshold')
print -dpdf ~/CS/BAC/results500closeFrom4500000using1eMinus4threshold.
plotRes('~/CS/BAC/dataForSim/simSec500_Nbacmix_500_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads.mat','~/CS/BAC/testSimSec500OriginalOrder')
print -dpdf ~/CS/BAC/results500closeFrom4500000usingOriginalOrder.pdf

for i=1:10
  subplot(2,5,i)
  plotRes('~/CS/BAC/dataForSim/simSec500_Nbacmix_500_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads.mat',['~/CS/BAC/testSimSec500_round',num2str(i)])
end
print -dpdf ~/CS/BAC/results500closeRepeatedPartitions

plotRes('~/CS/BAC/dataForSim/simSec500_Nbacmix_500_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads.mat','~/CS/BAC/testSimSec500_readLength100')
print -dpdf ~/CS/BAC/results500closeFrom4500000based_readlength100.pdf



%%%%%%%%%%%%%
% test original order
% run with original order
auxData = struct;
  
auxData.tmpFileName = ['testSimSec500OriginalOrder'];
auxData.keepOriginalOrderFlag = 1;
auxData.saveName = ['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/',auxData.tmpFileName];

auxData.createReadsFlag = 1;
auxData.randomizeSeedAfterCreateReadsFlag = 0;
auxData.brFlag = 0;
auxData.queueName = 'hour';
auxData.reads_data_fileName = '~/CS/BAC/dataForSim/simSec500_Nbacmix_500_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads.mat';
auxData.readLength = 50;
auxData.numProcessors = 100;
auxData.firstFlag = 0;
auxData.groupSize = 1000;
auxData.smallestSetOfCollected = 1000;
auxData.thresholdForCollectingBAC = 1e-3;
auxData.numBACtoConsider = 455055;
auxDataFileName = ['~/CS/tmp/structFor_',auxData.tmpFileName];
save(auxDataFileName,'auxData');
distributeBAC_general(auxDataFileName);

251227      251305      251537
252           252        252


%%%%%%%%%%%%%%555555
clear
%save ~/CS/BAC/dataForClose500 currInd ind_bac_in_mix correctWeight basicSeqNameDir basicSeqKey

load ~/CS/BAC/dataForClose500
readLength = 50;

%%%%% only 1000



%%%%%%%%%%%% 1000+1000
other = [1:455055];
other(ind_bac_in_mix) = [];
other = other(randperm(length(other)));
currData = [currInd,other(1:1000)];
[normalizedBac values] = prepareGroupOf1000DistributedSequenceFiles2(readLength,currData,basicSeqNameDir,basicSeqKey);  




currCorrect = correctWeight(currData)./sum(correctWeight(currData));
Yc = normalizedBac*currCorrect;
[currRes,firstFlag]=runOneGroupOf1000post(1,normalizedBac,Yc,1);

plot(currRes,currCorrect,'.')
title('1000+1000 randomly selected')
print -dpdf ~/CS/BAC/resClose500InTheir1000AndRandom1000


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% take only the 500 and add 1500
other = [1:455055];
other(ind_bac_in_mix) = [];
other = other(randperm(length(other)));
currData = [ind_bac_in_mix',other(1:1000)];
[normalizedBac values] = prepareGroupOf1000DistributedSequenceFiles2(readLength,currData,basicSeqNameDir,basicSeqKey);  

currCorrect = correctWeight(currData)./sum(correctWeight(currData));
Yc = normalizedBac*currCorrect;
[currRes,firstFlag]=runOneGroupOf1000post(1,normalizedBac,Yc,1);

plot(currRes,currCorrect,'.')
title('500close+1500 randomly selected')
print -dpdf ~/CS/BAC/resClose500AndRandom1500

find(currData==251537)


%%%%%%%%%%55
% take only those which are larger than 5*10^-3
other = [1:455055];
other(ind_bac_in_mix) = [];
other = other(randperm(length(other)));

a = find(correctWeight>5*10^-3);
currData = [a',other(1:1000)];

[normalizedBac values] = prepareGroupOf1000DistributedSequenceFiles2(readLength,currData,basicSeqNameDir,basicSeqKey);  

currCorrect = correctWeight(currData)./sum(correctWeight(currData));
Yc = normalizedBac*currCorrect;
[currRes,firstFlag]=runOneGroupOf1000post(1,normalizedBac,Yc,1);

plot(currRes,currCorrect,'.')
title('weight > 5*10^-3 +1000 randomly selected')
print -dpdf ~/CS/BAC/resCloseHigherThan5eMinus3AndRandom1000

% take only those which are larger than 2*10^-3
other = [1:455055];
other(ind_bac_in_mix) = [];
other = other(randperm(length(other)));

a = find(correctWeight>2*10^-3);
currData = [a',other(1:1000)];

[normalizedBac values] = prepareGroupOf1000DistributedSequenceFiles2(readLength,currData,basicSeqNameDir,basicSeqKey);  

currCorrect = correctWeight(currData)./sum(correctWeight(currData));
Yc = normalizedBac*currCorrect;
[currRes,firstFlag]=runOneGroupOf1000post(1,normalizedBac,Yc,1);

plot(currRes,currCorrect,'.')
title('weight > 2*10^-3 +1000 randomly selected')
print -dpdf ~/CS/BAC/resCloseHigherThan2eMinus3AndRandom1000

% take only those which are larger than 1.5*10^-3
other = [1:455055];
other(ind_bac_in_mix) = [];
other = other(randperm(length(other)));

a = find(correctWeight>1.5*10^-3);
currData = [a',other(1:1000)];

[normalizedBac values] = prepareGroupOf1000DistributedSequenceFiles2(readLength,currData,basicSeqNameDir,basicSeqKey);  

currCorrect = correctWeight(currData)./sum(correctWeight(currData));
Yc = normaylizedBac*currCorrect;
[currRes,firstFlag]=runOneGroupOf1000post(1,normalizedBac,Yc,1);

plot(currRes,currCorrect,'.')
title('weight > 1.5*10^-3 +1000 randomly selected')
print -dpdf ~/CS/BAC/resCloseHigherThan1point5eMinus3AndRandom1000


%%%%%%%%%%%%
% check for linear combinations
load ~/CS/BAC/dataForClose500
readLength = 50;

[normalizedBac values, ln] = prepareGroupOf1000DistributedSequenceFiles3(readLength,currInd,basicSeqNameDir,basicSeqKey);  

Cclose = normalizedBac'*normalizedBac;
pclose = symamd(Cclose);


tic;
[R,jb] = rref(normalizedBac(:,1:100));
toc;

junk = randperm(length(correctWeight));
junk = junk(1:1000);
[normalizedBac1 values,ln1] = prepareGroupOf1000DistributedSequenceFiles2(readLength,junk,basicSeqNameDir,basicSeqKey); 

Crandom = normalizedBac1'*normalizedBac1;
prandom = symamd(Crandom);

save ~/CS/BAC/linearInd500 Cclose pclose Crandom prandom normalizedBac1 normalizedBac
save ~/CS/BAC/linearInd500  normalizedBac1 normalizedBac ln ln1
%%%%%%%%

%%%%%%%%% solve twice: split into two groups for those which are close and others
imagesc(z)
o = linkage(1-z,'single');
c = cluster(o,'maxclust',400);

[junk,i1,i2] = unique(c);

take1 = currInd(i1);
currRes1 = solveSmall(readLength,take1,basicSeqNameDir,basicSeqKey,correctWeight);
a1 = find(currRes1>10^-3);

tmpC = c;
tmpC(i1) = [];
[junk,i1_2,i2] = unique(tmpC);
take2 = currInd(tmpC(i1_2));
currRes2 = solveSmall(readLength,take2,basicSeqNameDir,basicSeqKey,correctWeight);
a2 = find(currRes2>10^-3);


take = [take1(a1),take2(a2)];

currResTotal = solveSmall(readLength,take,basicSeqNameDir,basicSeqKey,correctWeight);

res = zeros(size(correctWeight));
res(take) = currResTotal;
plot(res,correctWeight,'.');


k=1;
for readLength=100:50:400
  currData = [currInd];
  [normalizedBac values] = prepareGroupOf1000DistributedSequenceFiles2(readLength,currData,basicSeqNameDir,basicSeqKey);  
  
  currCorrect = correctWeight(currData)./sum(correctWeight(currData));
  Yc = normalizedBac*currCorrect;
  [currRes,firstFlag]=runOneGroupOf1000post(1,normalizedBac,Yc,1);
  
  rnk = rank(full(normalizedBac));
  subplot(2,4,k)
  plot(currRes,currCorrect,'.')
  title(['readLength = ',num2str(readLength),' rank:',num2str(rnk)])
  k = k+1;
end
print -dpdf ~/CS/BAC/resClose500withDifferent_readLength


% summarize the whole thing 
two problems:

1. fn: freq in 1000 is 0 while should be positive
  add those which were not chosen but are close to those which were selected. maybe using the base
  
  candidate set: 
  maximize independence in groups of 1000 - in advance - randomly select and move those which are dependent from bucket to bucket
  threshold and the size of the gropu: set this as a function of number of reads that are explained by the set
  for each 1000 also save the values and distribution of their frequency
  
2. too many close cause singularity
  take 2 which are very similar and ask what should be done to
  differentiate between them
  show as a function oof read length
  
  locate region which differentite beteween them - lookat then and see the ratio.

3 linear combination or numerical issues
read length



why didn't we see it before

paired ends - prepare library


should we write?

deal with those do not appear in the database: 
   count how many reads are explained
   easy case:those which appear are far
   check the histogram of read for each bacteria: are these consistent? 


simulations: change code, test: several 1000 close 
    
%%%%%%%%%%%%%%%%%% try differnet 

options = optimset('Display','iter','TolFun',1e-8);
newX = lsqlin(normalizedBac,Yc,[],[],ones(1,1000),1);

newX = lsqnonneg(normalizedBac,Yc,10^-10);


options = optimset('Display','iter','TolFun',1e-12);
newX = lsqnonneg(normalizedBac,Yc,options);

[x{i},firstFlag]=runOneGroupOf1000postTmp(1,normalizedBac,Yc,1);

%%%%%%%%%%% add 10000 to these 500



currCorrect = correctWeight(currInd)./sum(correctWeight(currInd)); 
q = ciel(normalizedBac);
a = find(sum(q,2)>500);
el = 1:size(normalizedBac,1);
el(a) = [];
Yca = Yc;
Yca(a) = [];
[currRes,firstFlag]=runOneGroupOf1000post(1,normalizedBac(el,:),Yca,1);
plot(currRes,currCorrect,'.')


currResJust500 = solveSmall(readLength,ind_bac_in_mix,basicSeqNameDir,basicSeqKey,correctWeight);
res = zeros(size(correctWeight));
res(ind_bac_in_mix) = currResJust500;
plot(res,correctWeight,'.');


[normalizedBac values] = prepareGroupOf1000DistributedSequenceFiles2(readLength,ind_bac_in_mix,basicSeqNameDir,basicSeqKey);  

A = normalizedBac;
A(end+1,:) = 1/500;

currCorrect = correctWeight(ind_bac_in_mix)./sum(correctWeight(ind_bac_in_mix));
YcA = [normalizedBac*currCorrect;1];

currResJust500 = YcA\A;
res = zeros(size(correctWeight));
res(ind_bac_in_mix) = currResJust500;
plot(res(ind_bac_in_mix),correctWeight(ind_bac_in_mix),'.');


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example of 1000 close 
clear 
load ~/CS/BAC/dataForClose500
basicSeqKey = '/homes/csfaculty/shental/CS/BAC/dat450000/key450000';
basicSeqNameDir = '/homes/csfaculty/shental/CS/BAC/dat450000/';


readLength = 50;

% solve using independent
uniqueReads = 1;
uniqueReads_length = 1;
[normalizedBac values independentInd{1},dependentInd{1}] = prepareGroupOf1000DistributedSequenceFiles4(readLength,currInd,basicSeqNameDir,basicSeqKey);

[junk,i1,i2] = intersect(ind_bac_in_mix,independentInd{1});
currCorrect = zeros(size(normalizedBac,2),1);
currCorrect(i2) = correctWeight(ind_bac_in_mix(i1))';
Yca = normalizedBac*currCorrect;
[currRes,firstFlag]=runOneGroupOf1000post(1,normalizedBac,Yca,0);

tmpInd{1} = currInd;
[x_extra,ind_extra] = dealWithDependent5(dependentInd,5,[],tmpInd,basicSeqNameDir,basicSeqKey,readLength,uniqueReads,uniqueReads_length);


clf
plot(currRes,correctWeight(independentInd{1})','.')
hold on
plot(x_extra,correctWeight(ind_extra)','ro')

% solve the whole 1000
[normalizedBac2 values2 ] = prepareGroupOf1000DistributedSequenceFiles2(readLength,currInd,basicSeqNameDir,basicSeqKey);
currCorrect = correctWeight(currInd)./sum(correctWeight(currInd));
Yca2 = normalizedBac2*currCorrect;;
[currRes2,firstFlag]=runOneGroupOf1000post(1,normalizedBac2,Yca2,1);
figure;
plot(currRes2,currCorrect,'.')



%%%%%%%%%%%%55
clear
load ~/CS/BAC/dataForClose500
basicSeqKey = '/homes/csfaculty/shental/CS/BAC/dat450000/key450000';
basicSeqNameDir = '/homes/csfaculty/shental/CS/BAC/dat450000/';

correctWeight1 = correctWeight;
w1 = correctWeight(ind_bac_in_mix);
clear correctWeight w


Nbac_in_mixture = 500;
correctWeight = correctWeight1;
ind_bac_in_mix = ind_bac_in_mix;
allBAC = currInd;
outseed = [];
Nreads = 10^6;
readlen = 50;
bac_dist_flag = 1;
npower = 0.5  ;         
w = w1;

save /homes/csfaculty/shental/CS/BAC/dataForSim/sim500CloseDifferentCases_Nbacmix_1000_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads  Nbac_in_mixture correctWeight    outseed    Nreads   ind_bac_in_mix   readlen      bac_dist_flag    npower    w allBAC
     
clear
strinfilename = '/homes/csfaculty/shental/CS/BAC/dataForSim/sim500CloseDifferentCases_Nbacmix_1000_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads';
                 
auxData = struct;
auxData.inifiniteNumberOfReadsFlag = 0;
auxData.reads_data_fileName = strinfilename;
auxData.minimalNumberOfDependent = 1;
auxData.tmpFileName = 'test1000CloseFinite';
auxData.currDir = ['/homes/csfaculty/shental/CS/BAC/tmpRuns/',auxData.tmpFileName];
unix(['mkdir ',auxData.currDir]);
unix(['cd ',auxData.currDir]); 
  
auxData.readLength = 50;
auxData.inifiniteNumberOfReadsFlag = 0;
auxData.keepOriginalOrderFlag = 1;
auxData.saveName = ['~/CS/BAC/',auxData.tmpFileName];
auxData.numProcessors = 1;
auxData.groupSize = 999;
auxData.smallestSetOfCollected = 999;
auxData.numBACtoConsider = 1000;
auxData.thresholdForCollectingBAC = 1e-3;
auxData.firstFlag = 1;
auxData.createReadsFlag = 1;
auxDataFileName = ['~/CS/tmp/structFor_',auxData.tmpFileName];
save(auxDataFileName,'auxData')
distributeBAC_generalTmp(auxDataFileName)

load ~/CS/BAC/dataForSim/sim500CloseDifferentCases_Nbacmix_1000_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads
load ~/CS/BAC/test1000CloseInFinite

figure(1)
for i=1:length(found)
  if ~isempty(found{i})
    
    plot(found{i}(ind_bac_in_mix),correctWeight(ind_bac_in_mix)','.')
  end
end

load ~/CS/BAC/dataForSim/sim500CloseDifferentCases_Nbacmix_1000_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads
load ~/CS/BAC/test1000CloseInFinite1

figure(2)
for i=1:length(found)
  if ~isempty(found{i})
    
    plot(found{i}(ind_bac_in_mix),correctWeight(ind_bac_in_mix)','.')
  end
end

load ~/CS/BAC/dataForSim/sim500CloseDifferentCases_Nbacmix_1000_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads
load ~/CS/BAC/test1000CloseFinite

figure(3)
for i=1:length(found)
  if ~isempty(found{i})
    
    plot(found{i}(ind_bac_in_mix),correctWeight(ind_bac_in_mix)','.')
  end
end

load ~/CS/BAC/dataForSim/sim500CloseDifferentCases_Nbacmix_1000_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads
load ~/CS/BAC/test1000CloseFinite1

figure(4)
for i=1:length(found)
  if ~isempty(found{i})
    
    plot(found{i}(ind_bac_in_mix),correctWeight(ind_bac_in_mix)','.')
  end
end


plotRes('~/CS/BAC/dataForSim/sim500CloseDifferentCases_Nbacmix_1000_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads','~/CS/BAC/test1000CloseFinite')

do we change solve again?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555
% change the 

clear
i = 9;
load ~/CS/BAC/dataForSim/simSec500_Nbacmix_500_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads
load(['~/CS/BAC/testSimSec500_round',num2str(i)]);
r = found{end};

a = find(r==0 & correctWeight'>0.01);
length(store_kp)

for j=1:length(store_kp)
  if isempty(find(store_kp{j}==a))
    disp(['did not find it in: ',num2str(j)])
    b = find(store_kp{j-1}==a);
    store_X{j-1}(b)
  end
end



