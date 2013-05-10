%%%%%%%%%%
clear
auxData = struct;
auxData.outFileDirName = 'canonical10e3';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [1]; % power law 1
auxData.readLength = [100];
auxData.numIter = 10;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [0 1];
auxData.correctMeasurementForReadErrorsFlagInput = [0 1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [50];
auxData.numProcessors{2}([50 100 200 400]) = [50 100 75 150];
auxData.doL1Flag = 1;
runGeneralSplitCorrection(auxData);

clear
auxData = struct;
auxData.outFileDirName = 'canonical10e3_onlyNoiseCorrected';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [1]; % power law 1
auxData.readLength = [100];
auxData.numIter = 10;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [50];
auxData.numProcessors{2}([50 100 200 400]) = [50 100 75 150];
auxData.doL1Flag = 1;
runGeneralSplitCorrection(auxData);

clear
auxData = struct;
auxData.outFileDirName = 'canonical10e3_full';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [1]; % power law 1
auxData.readLength = [100];
auxData.numIter = 50;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [0 1];
auxData.correctMeasurementForReadErrorsFlagInput = [0 1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10];
auxData.numProcessors{2}([50 100 200 400]) = [50 50 75 150];
auxData.doL1Flag = 1;
runGeneralSplitCorrection(auxData);

clear
auxData = struct;
auxData.outFileDirName = 'canonical10e3_full2';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [1]; % power law 1
auxData.readLength = [100];
auxData.numIter = 50;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [0 1];
auxData.correctMeasurementForReadErrorsFlagInput = [0 1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10];
auxData.numProcessors{2}([50 100 200 400]) = [50 50 75 150];
auxData.doL1Flag = 1;
runGeneralSplitCorrection(auxData);


% %%%%%%%%%
% read length
clear
auxData = struct;
auxData.outFileDirName = 'VariableReadLength';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [1]; % power law 1
auxData.readLength = [35 100 200];
auxData.numIter = 10;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10];
auxData.numProcessors{2}([35 50 100 200 400]) = [50 50 50 75 150];
auxData.doL1Flag = 1;
runGeneralSplitCorrection(auxData);


% read length
clear
auxData = struct;
auxData.outFileDirName = 'VariableReadLength_200';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [1]; % power law 1
auxData.readLength = [200];
auxData.numIter = 10;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10];
auxData.numProcessors{2}([35 50 100 200 400]) = [50 50 50 100 150];
auxData.doL1Flag = 1;
runGeneralSplitCorrection(auxData);

% read length full
clear
auxData = struct;
auxData.outFileDirName = 'VariableReadLength_full';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [1]; % power law 1
auxData.readLength = [200 35];
auxData.numIter = 50;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10];
auxData.numProcessors{2}([35 50 100 200 400]) = [50 50 50 100 150];
auxData.doL1Flag = 1;
runGeneralSplitCorrection(auxData);

% read length full2
clear
auxData = struct;
auxData.outFileDirName = 'VariableReadLength_full2';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [1]; % power law 1
auxData.readLength = [200 35];
auxData.numIter = 50;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10];
auxData.numProcessors{2}([35 50 100 200 400]) = [50 50 50 100 150];
auxData.doL1Flag = 1;
runGeneralSplitCorrection(auxData);



% number of bacteria
clear
auxData = struct;
auxData.outFileDirName = 'VariableNumberOfBacteria';
auxData.bdf = 0;
auxData.nb = [10 100 300 500];
auxData.pwr = [1]; % power law 1
auxData.readLength = [100];
auxData.numIter = 10;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10];
auxData.numProcessors{2}([50 100 200 400]) = [50 50 75 150];
auxData.doL1Flag = 1;
runGeneralSplitCorrection(auxData);


clear
auxData = struct;
auxData.outFileDirName = 'VariableNumberOfBacteria_2';
auxData.bdf = 0;
auxData.nb = [200 400 700];
auxData.pwr = [1]; % power law 1
auxData.readLength = [100];
auxData.numIter = 10;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10];
auxData.numProcessors{2}([50 100 200 400]) = [50 50 75 150];
auxData.doL1Flag = 1;
runGeneralSplitCorrection(auxData);

clear
auxData = struct;
auxData.outFileDirName = 'VariableNumberOfBacteria_full';
auxData.bdf = 0;
auxData.nb = [10 100:200:1000];
auxData.pwr = [1]; % power law 1
auxData.readLength = [100];
auxData.numIter = 50;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10];
auxData.numProcessors{2}([50 100 200 400]) = [50 50 75 150];
auxData.doL1Flag = 1;
runGeneralSplitCorrection(auxData);

clear
auxData = struct;
auxData.outFileDirName = 'VariableNumberOfBacteria_full2';
auxData.bdf = 0;
auxData.nb = [10 100:200:1000];
auxData.pwr = [1]; % power law 1
auxData.readLength = [100];
auxData.numIter = 50;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10];
auxData.numProcessors{2}([50 100 200 400]) = [50 50 75 150];
auxData.doL1Flag = 1;
runGeneralSplitCorrection(auxData);



%%%%%%%%%%%%%%% end number of bacteria


%%%%%%%%%%55
% number of reads 0- should be done again as for less than 10^5 reads - there was no noise added
clear
auxData = struct;
auxData.outFileDirName = 'VariableNumberOfReads';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [1]; % power law 1
auxData.readLength = [100];
auxData.numIter = 10;
auxData.nr = [10^4 5*10^5];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10];
auxData.numProcessors{2}([50 100 200 400]) = [50 100 75 150];
auxData.doL1Flag = 1;
runGeneralSplitCorrection(auxData);

%%%%%%%%%%%%%%%%%%%%%%%
clear
auxData = struct;
auxData.outFileDirName = 'VariableNumberOfReads_10e4';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [1]; % power law 1
auxData.readLength = [100];
auxData.numIter = 10;
auxData.nr = [10^4];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10];
auxData.numProcessors{2}([50 100 200 400]) = [50 100 75 150];
auxData.doL1Flag = 1;
runGeneralSplitCorrection(auxData);


clear
auxData = struct;
auxData.outFileDirName = 'VariableNumberOfReads_full';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [1]; % power law 1
auxData.readLength = [100];
auxData.numIter = 50;
auxData.nr = [10^4 10^5 5*10^5];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10];
auxData.numProcessors{2}([50 100 200 400]) = [50 100 75 150];
auxData.doL1Flag = 1;
runGeneralSplitCorrection(auxData);

