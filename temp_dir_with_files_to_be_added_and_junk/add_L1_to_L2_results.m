clear
auxData = struct;
auxData.outFileDirName = 'numberOfBacteria';
auxData.bdf = 0;
auxData.nb = [10 100 500 1000];
auxData.pwr = [0];
auxData.readLength = [100];
auxData.numIter = 50;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [0 1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors = [10 35];

% this can probably be done only with no_noise case

clear
basicName = 'canonicalFull_bac_dist_flag_0_Nbac_in_mixture_1000_npower_0_readlen_100_numIter_9_Nreads_1000000_Noise_0_Correction_0';

auxDataFileName = ['/home/csfaculty/shental/CS/tmp/canonicalFull/structFor_',basicName];
resFile = ['/home/csfaculty/shental/CS/BAC/canonicalFull/',basicName];

check_L1_finalStage_uniformCoverageSim(basicName,auxDataFileName,resFile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create new run to compare L1 and L2

clear
auxData = struct;
auxData.outFileDirName = 'forCompareL1L2';
auxData.bdf = 0;
auxData.nb = [50];
auxData.pwr = [0];
auxData.readLength = [100];
auxData.numIter = 2;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [0];
auxData.correctMeasurementForReadErrorsFlagInput = [0];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [50];
auxData.numProcessors{2}([50 100 200 400]) = [25 50 75 150];
runGeneralSplitCorrection(auxData);


% move the results in br
moverResults(auxData.outFileDirName);

% move to place in sol
extractResults(auxData.outFileDirName)

% check the results
clear
basicName = 'forCompareL1L2_bac_dist_flag_0_Nbac_in_mixture_50_npower_0_readlen_100_numIter_1_Nreads_1000000_Noise_0_Correction_0';

auxDataFileName = ['/home/csfaculty/shental/CS/tmp/forCompareL1L2/structFor_',basicName];
resFile = ['/home/csfaculty/shental/CS/BAC/forCompareL1L2/',basicName];

check_L1_finalStage_uniformCoverageSim(basicName,auxDataFileName,resFile)

