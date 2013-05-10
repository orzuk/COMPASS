clear
resultsDir = ['~/CS/BAC/12Samples/Solexa/results'];
userDir = getuserdir;

auxDataTest.datasetName = 'primers750_primer_tail';
auxDataTest.numBACtoConsider = 212040;
auxDataTest.dataSetPrimers750Flag = 1;
auxDataTest.readLength = 90;
auxDataTest.basicDirNameForFile = ['window_',num2str(auxDataTest.readLength),'_withNor_MergeForwardAndReverse'];

% load the reads from unnormalized 
readDataFileName = [userDir,'/CS/BAC/12Samples/Solexa/data/window_90_noNor_MergeForwardAndReverse/data_window_90_noNor_MergeForwardAndReverse_S1'];
outputName = 'test_window_90_withNor_MergeForwardAndReverse_normReadsf_S1';
profileName = '~/CS/BAC/12Samples/Solexa/data/GlobalReadProfile/GlobReadProfileS1';
extractRes_weightedMatrix(outputName,readDataFileName,profileName,resultsDir,auxDataTest,1)


