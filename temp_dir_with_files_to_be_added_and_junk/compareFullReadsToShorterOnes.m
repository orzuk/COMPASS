% create reads
% add noise - regular
clear
auxData = struct;
auxData.outFileDirName = 'test400VscanonicalFull';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [0 0.5];
auxData.readLength = [100];
auxData.numIter = 10;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [0 1];
auxData.correctMeasurementForReadErrorsFlagInput = [0 1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors = [10 35];
runGeneralSplitCorrection(auxData);

% add random noise along the 350


% run with noise, without noise