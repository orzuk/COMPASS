if 1==2

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
     
% finite
clear
strinfilename = '/homes/csfaculty/shental/CS/BAC/dataForSim/sim500CloseDifferentCases_Nbacmix_1000_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads';
                 
auxData = struct;
auxData.inifiniteNumberOfReadsFlag = 0;
auxData.reads_data_fileName = strinfilename;
auxData.minimalNumberOfDependent = 1;
auxData.tmpFileName = 'test1000CloseFiniteNormalize';
auxData.normalizeByTotalNumberOfReadsFlag = 0;
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


% infinite
clear
strinfilename = '/homes/csfaculty/shental/CS/BAC/dataForSim/sim500CloseDifferentCases_Nbacmix_1000_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads';
                 
auxData = struct;
auxData.inifiniteNumberOfReadsFlag = 1;
auxData.reads_data_fileName = strinfilename;
auxData.minimalNumberOfDependent = 1;
auxData.tmpFileName = 'test1000CloseInFinite';
auxData.currDir = ['/homes/csfaculty/shental/CS/BAC/tmpRuns/',auxData.tmpFileName];
%unix(['mkdir ',auxData.currDir]);
%unix(['cd ',auxData.currDir]); 
  
auxData.readLength = 50;
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



end

% check the effect of finite numberof reads
figure(2)
subplot(2,2,1)
plotClose('~/CS/BAC/test1000CloseInFinite',1)
subplot(2,2,2)
plotClose('~/CS/BAC/test1000CloseInFinite1',2)
subplot(2,2,3)
plotClose('~/CS/BAC/test1000CloseFinite',3)
subplot(2,2,4)
plotClose('~/CS/BAC/test1000CloseFinite1',4)


% conmpare finite and infinite


% check the effect of normalizing by the local number of reads
figure(1)
subplot(2,2,1)
plotClose('~/CS/BAC/test1000CloseFiniteDoNotNormalize.mat',1)
subplot(2,2,2)
plotClose('~/CS/BAC/test1000CloseFiniteDoNotNormalize1.mat',2)
subplot(2,2,3)
plotClose('~/CS/BAC/test1000CloseFiniteNormalize.mat',3)
subplot(2,2,4)
plotClose('~/CS/BAC/test1000CloseFiniteNormalize1.mat',4)




