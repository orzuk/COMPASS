% run the data - before normalization
clear
basicName = {'S1','S2','S3','S4','M1','M2','M3','M4','S7','S10','O7','O10'};

%%%%%


%%%%%%%%


auxDataTest.datasetName = 'primers750_primer_tail';
auxDataTest.numBACtoConsider = 212040;
auxDataTest.dataSetPrimers750Flag = 1;

auxDataTest.readLength = 100;

auxDataTest.currNumProcessors = 20;
auxDataTest.basicDirNameForFile = ['regular_',num2str(auxDataTest.readLength),'_withNor'];

clear fileName inputName outputName
for i=1:length(basicName)
  fileName{i} = ['data_',auxDataTest.basicDirNameForFile,'_',basicName{i}];
  outputName{i} = ['test_',auxDataTest.basicDirNameForFile,'_',basicName{i}];
end

% prepare the files

if 1==2
  % multiply the reads by 10000
  for i=1:length(basicName)
    dat = load([userDir,'/CS/BAC/12Samples/Solexa/data/',auxDataTest.basicDirNameForFile,'/NormReads',basicName{i}]);
    datUnique = load([userDir,'/CS/BAC/12Samples/Solexa/data/opt3_th50_mar25_2012/readsSolexa_with_equalization_opt3_th50_mar25_2012_100_',basicName{i}]);
    
    clear uniqueReads uniquereads_length
    uniqueReads = datUnique.uniqueReads;
    uniqueReads_length = ceil(dat.normReadsf*10^6);
    
    length(find(uniqueReads_length==1))
    
    save([userDir,'/CS/BAC/12Samples/Solexa/data/',auxDataTest.basicDirNameForFile,'/',fileName{i}],'uniqueReads','uniqueReads_length')
  end
end

% copy to br

% prepare the files
for i=1%:length(outputName)
  testSimFromReadsLikeRealData_generic2(outputName{i},['12Samples/Solexa/data/',auxDataTest.basicDirNameForFile,'/',fileName{i}],auxDataTest)
end

%%%%%5
% run the files

for i=1%:length(outputName)
  unix(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_',outputName{i},'/run_sol_',outputName{i},'_noCorrection_',num2str(auxDataTest.readLength)])
end

% collect results
unix(['mkdir /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/12Samples/Solexa/results/',auxDataTest.basicDirNameForFile])
for i=1%:length(outputName)
  unix(['cp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/sol_',outputName{i},'/sol_',outputName{i},'_noCorrection_',num2str(auxDataTest.readLength),'/sol_',outputName{i},'_noCorrection_',num2str(auxDataTest.readLength),'.mat',' /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/12Samples/Solexa/results/',auxDataTest.basicDirNameForFile])
end


unix(['cd /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/12Samples/Solexa/results/',';tar cf results_', ...
      auxDataTest.basicDirNameForFile,'.tar ',auxDataTest.basicDirNameForFile])

% copy locally
['scp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/12Samples/Solexa/results/results_',auxDataTest.basicDirNameForFile,'.tar shental@sol.cslab.openu.ac.il:~/CS/BAC/12Samples/Solexa/results/']
% run locally

resultsDir = ['~/CS/BAC/12Samples/Solexa/results'];
for i=1:length(outputName)
  readDataFileName = [userDir,'/CS/BAC/12Samples/Solexa/data/',auxDataTest.basicDirNameForFile,'/',fileName{i}]
  extractResRealData_generic(outputName{i},readDataFileName,resultsDir,auxDataTest)
  drawnow
  %pause
end


for i=1%:length(outputName)
  plotResRealData_generic(resultsDir,outputName{i},auxDataTest)
  drawnow
  pause
end


