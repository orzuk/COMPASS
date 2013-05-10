clear
auxData = struct;
auxData.outFileDirName = 'canonicalFormWithL1InLastStage';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [0 0.5];
auxData.readLength = [100];
auxData.numIter = 1;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [0 1];
auxData.correctMeasurementForReadErrorsFlagInput = [0 1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [50];
auxData.numProcessors{2}([50 100 200 400]) = [100 100 75 150];
auxData.doL1Flag = 1;
runGeneralSplitCorrection(auxData);

%numProcessors = 100;
%dirName = 'canonicalFormWithL1InLastStage';
%changeNumberOfProcessors(dirName,numProcessors)

dir = [];
solveL1After(auxData.outFileDirName,0)

%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
auxData = struct;
auxData.outFileDirName = 'canonicalFormWithL1';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [0.5];
auxData.readLength = [100];
auxData.numIter = 10;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [0 1];
auxData.correctMeasurementForReadErrorsFlagInput = [0 1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [20];
auxData.numProcessors{2}([50 100 200 400]) = [50 50 75 150];
auxData.doL1Flag = 1;
runGeneralSplitCorrection(auxData);

dir = [];
solveL1After(auxData.outFileDirName,4,1:29)


clear
auxData = struct;
auxData.outFileDirName = 'canonicalFormWithL1_HighCoverage';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [0.5];
auxData.readLength = [100];
auxData.numIter = 5;
auxData.nr = [10^7];
auxData.addNoiseFlagInput = [0 1];
auxData.correctMeasurementForReadErrorsFlagInput = [0 1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;

auxData.numProcessors{1} = [50];
auxData.numProcessors{2}([50 100 200 400]) = [50 100 75 150];
auxData.doL1Flag = 1;
runGeneralSplitCorrection(auxData);

% move tmp and dataForSim to sol

% run in sol then transfer the reads to br
%userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
userDir = getuserdir;
doLocalCreateReadsForSpecificCasesGeneric(auxData,10^6,100,userDir)

% move to br
numProcessors = 100;
changeNumberOfProcessors(auxData.outFileDirName,numProcessors)

% the corrected did not work - took too long to correct even 1 
solveL1After(auxData.outFileDirName,4,1:9)


%%%%%%%%%%%%5555
% 
clear
auxData = struct;
auxData.outFileDirName = 'tsl_local';
auxData.bdf = 0;
auxData.nb = [10];
auxData.pwr = [0.5];
auxData.readLength = [100];
auxData.numIter = 1;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [0];
auxData.correctMeasurementForReadErrorsFlagInput = [0];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [50];
auxData.numProcessors{2}([50 100 200 400]) = [50 50 75 150];
auxData.doL1Flag = 1;
runGeneralSplitCorrection(auxData);
