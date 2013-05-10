clear
basicName = {'S7','S10','O7','O10'};

auxDataTest.datasetName = 'data750with4mm';
auxDataTest.numBACtoConsider = 299783;
auxDataTest.dataSetPrimers750Flag = 1;

auxDataTest.readLength = 90;

auxDataTest.currNumProcessors = 30;
auxDataTest.basicDirNameForFile = ['window_',num2str(auxDataTest.readLength),'_withNor_MergeForwardAndReverse_data750withMM'];


clear fileName inputName outputName
for i=1:length(basicName)
  %fileName{i} = ['data_',auxDataTest.basicDirNameForFile,'_normReadsf_',basicName{i}];
  fileName{i} = ['data_window_90_withNor_MergeForwardAndReverse_normReadsf_',basicName{i}];
  outputName{i} = ['test_',auxDataTest.basicDirNameForFile,'_normReadsf_',basicName{i}];
end

% prepare the files
for i=1:length(outputName)
  testSimFromReadsLikeRealData_generic2(outputName{i},['12Samples/Solexa/data/window_90_withNor_MergeForwardAndReverse/',fileName{i}],auxDataTest)
end

%%%%%5
% run the files

for i=1:length(outputName)
  unix(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_',outputName{i},'/run_sol_',outputName{i},'_noCorrection_',num2str(auxDataTest.readLength)])
end

% collect results
unix(['mkdir /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/12Samples/Solexa/results/',auxDataTest.basicDirNameForFile])
for i=1:length(outputName)
  unix(['cp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/sol_',outputName{i},'/sol_',outputName{i},'_noCorrection_',num2str(auxDataTest.readLength),'/sol_',outputName{i},'_noCorrection_',num2str(auxDataTest.readLength),'.mat',' /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/12Samples/Solexa/results/',auxDataTest.basicDirNameForFile])
end
unix(['cd /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/12Samples/Solexa/results/',';tar cf results_', ...
      auxDataTest.basicDirNameForFile,'.tar ',auxDataTest.basicDirNameForFile])

% copy locally
['scp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/12Samples/Solexa/results/results_',auxDataTest.basicDirNameForFile,'.tar shental@sol.cslab.openu.ac.il:~/CS/BAC/12Samples/Solexa/results/']
% run locally



resultsDir = ['~/CS/BAC/12Samples/Solexa/results'];
userDir = getuserdir;
for i=1:length(outputName)
  readDataFileName = [userDir,'/CS/BAC/12Samples/Solexa/data/window_90_withNor_MergeForwardAndReverse/',fileName{i}]
  extractResRealData_generic(outputName{i},readDataFileName,resultsDir,auxDataTest)
  drawnow
  %pause
end

%%%%%%%%%%%5555
% load the reads from unnormalized 
filterFlag = 1;
for i=1:length(outputName)
  readDataFileName_noNor = [userDir,'/CS/BAC/12Samples/Solexa/data/window_90_noNor_MergeForwardAndReverse/data_window_90_noNor_MergeForwardAndReverse_',basicName{i}];
  profileName = [userDir,'/CS/BAC/12Samples/Solexa/data/GlobalReadProfile/Noise',basicName{i}];
  extractRes_weightedMatrix(outputName{i},readDataFileName_noNor,profileName,resultsDir,auxDataTest,filterFlag)
end
 
 




