function createRun_simClosePeter(mixtureName,Nreads,currNumProcessors)
%keyboard
if  ~isempty(findstr(getenv('HOSTNAME'),'openu'))
  userDir = getuserdir;
else
  userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
end


readLength = 100;
addNoiseFlag = 1;

auxDataTest = struct;
auxDataTest.datasetName = 'fullDatabse';
auxDataTest.numBACtoConsider = 410849;
auxDataTest.dataSetPrimers750Flag = 0;
auxDataTest.readLength = 100;
auxDataTest.currNumProcessors = currNumProcessors; % change for the with correction case
auxDataTest.correctMeasurementForReadErrorsFlag = 1; % run with correction
dataName = [mixtureName,'_reads_readlen_',num2str(readLength),'_noise_',num2str(addNoiseFlag),'_Nreads_',num2str(Nreads)];

basicDirNameForFile = [dataName,'_correction_',num2str(auxDataTest.correctMeasurementForReadErrorsFlag)]
outputName = ['test_',basicDirNameForFile];

testLikeRealData_genericChooseCorrection(outputName,['toyExample/',mixtureName,'/reads/',dataName],auxDataTest)

unix(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_',outputName,'/run_sol_',outputName,'_withCorrection_readlen_',num2str(auxDataTest.readLength)])

% copy the results to one directoryx


