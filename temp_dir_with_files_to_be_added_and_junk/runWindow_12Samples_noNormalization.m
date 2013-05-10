% run the data - before normalization
clear
basicName = {'S1','S2','S3','S4','M1','M2','M3','M4','S7','S10','O7','O10'};

auxDataTest.datasetName = 'primers750_primer_tail';
auxDataTest.numBACtoConsider = 212040;
auxDataTest.dataSetPrimers750Flag = 1;

auxDataTest.readLengthPre = 100;
auxDataTest.readLengthPost = 90;
auxDataTest.readLength = auxDataTest.readLengthPost;

auxDataTest.currNumProcessors = 20;
auxDataTest.basicDirNameForFile = ['window_',num2str(auxDataTest.readLengthPost),'_noNor'];

clear fileName inputName outputName
for i=1:length(basicName)
  fileName{i} = ['data_',auxDataTest.basicDirNameForFile,'_',basicName{i}];
  outputName{i} = ['test_',auxDataTest.basicDirNameForFile,'_',basicName{i}];
end


% prepare the files

for i=1:length(outputName)
  testSimFromReadsLikeRealData_generic2(outputName{i},['12Samples/Solexa/data/',auxDataTest.basicDirNameForFile,'/',fileName{i}],auxDataTest)
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
unix(['scp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/12Samples/Solexa/results/results_',auxDataTest.basicDirNameForFile,'.tar shental@sol.cslab.openu.ac.il:~/CS/BAC/12Samples/Solexa/results/'])
% run locally


for i=1:length(outputName)
  extractResSim_generic(outputName{i},inputName{i},fileName{i},auxDataTest)
  drawnow
  %pause
end

for i=[1:10 12:15]%length(outputName)
  plotResSim_generic(outputName{i},inputName{i},fileName{i},auxDataTest)
  drawnow
  pause
end


