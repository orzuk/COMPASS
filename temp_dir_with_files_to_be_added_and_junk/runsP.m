clear
auxData = struct;
auxData.outFileDirName = 'canonicalFull';
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

changeNumberOfProcessors(dirName,numProcessors)
% was done for 0 not 0.5
% first 4 started 14:30

clear
auxData = struct;
auxData.outFileDirName = 'canonicalFull_add40';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [0 0.5];
auxData.readLength = [100];
auxData.numIter = 40;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [0 1];
auxData.correctMeasurementForReadErrorsFlagInput = [0 1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors = [10 35];
runGeneralSplitCorrection(auxData);

changeNumberOfProcessors('canonicalFull_add40',50)
% run list1 notice to run only the 0 and not 0.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear
auxData = struct;
auxData.outFileDirName = 'canonical750';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [0 0.5];
auxData.readLength = [100];
auxData.numIter = 10;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [0 1];
auxData.correctMeasurementForReadErrorsFlagInput = [0 1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 1;
auxData.numProcessors = [10 35];
runGeneralSplitCorrection(auxData);

changeNumberOfProcessors(dirName,numProcessors)
changeList([userDir,/CS/BAC/dataForSim/numberOfBacteria/])
clear
auxData = struct;
auxData.outFileDirName = 'canonical750_add40';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [0 0.5];
auxData.readLength = [100];
auxData.numIter = 40;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [0 1];
auxData.correctMeasurementForReadErrorsFlagInput = [0 1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 1;
auxData.numProcessors = [10 35];
runGeneralSplitCorrection(auxData);

changeNumberOfProcessors(dirName,numProcessors)

%%
clear
auxData = struct;
auxData.outFileDirName = 'numberOfReads';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [0];
auxData.readLength = [100];
auxData.numIter = 50;
auxData.nr = [10^4 5*10^4 10^5 5*10^5 10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [0 1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors = [10 35];
runGeneralSplitCorrection(auxData);

changeNumberOfProcessors(dirName,numProcessors)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
runGeneralSplitCorrection(auxData);

changeNumberOfProcessors(dirName,numProcessors)



changeNumberOfProcessors('canonical750',50);changeNumberOfProcessors('canonical750_add40',50);changeNumberOfProcessors('numberOfReads',50);changeNumberOfProcessors('numberOfBacteria',50);

[junk,val] = unix(['bjobs -u orzuk | grep hour | grep RUN | wc ']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
auxData = struct;
auxData.outFileDirName = 'numberOfBacteriaRerun10';
auxData.bdf = 0;
auxData.nb = [10];
auxData.pwr = [0];
auxData.readLength = [100];
auxData.numIter = 100;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors = [10 35];
runGeneralSplitCorrection(auxData);
changeNumberOfProcessors(auxData.outFileDirName,50)

%%%%%%%%%%%%%%%%%%%%%%%%55
clear
auxData = struct;
auxData.outFileDirName = 'testReadLength';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [0];
auxData.readLength = [400];%50 100 150 200 300 400
auxData.numIter = 1;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors = [10 100];
runGeneralSplitCorrection(auxData);

userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
doLocalCreateReadsForSpecificCases(auxData,userDir)

%%%%%%%%%%%%%%%%%%%%%%%%55

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
auxData = struct;
auxData.outFileDirName = 'varyReadLength';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [0];
auxData.readLength = [50 100 200 400];
auxData.numIter = 10;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10];
auxData.numProcessors{2}([50 100 200 400]) = [25 50 75 150];
runGeneralSplitCorrection(auxData);

userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
doLocalCreateReadsForSpecificCases(auxData,userDir)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
auxData = struct;
auxData.outFileDirName = 'varyReadLength2';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [0];
auxData.readLength = [50 100 200 400];
auxData.numIter = 10;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10];
auxData.numProcessors{2}([50 100 200 400]) = [25 50 75 150];
runGeneralSplitCorrection(auxData);

userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
doLocalCreateReadsForSpecificCases(auxData,userDir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
auxData = struct;
auxData.outFileDirName = 'varyReadLength3';
auxData.bdf = 0;
auxData.nb = [1000];
auxData.pwr = [0];
auxData.readLength = [50 100 200 400];
auxData.numIter = 30;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10];
auxData.numProcessors{2}([50 100 200 400]) = [25 50 75 150];
runGeneralSplitCorrection(auxData);

userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
doLocalCreateReadsForSpecificCases(auxData,userDir)

