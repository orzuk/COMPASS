function runDataPeter_no_eq_with_correction(basicDirName,readsDataFile)

disp([readsDataFile])
readLength = 101;
currNumProcessors = 200;

numBACtoConsider = 410849;
sampleName{1} = basicDirName;

for i=1
  
  clear auxData
  
  %userDir = getuserdir;
  userDir =        '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
  addpath([userDir,'/CS/mFilesBAC/']);
  addpath(genpath([userDir,'/CS/BAC/cvx/']))
  
  dirName = ['sol_',basicDirName,'_',sampleName{i}];
  
  
  unix(['mkdir ',userDir,'/CS/tmp/',dirName,...
        ';mkdir ',userDir,'/CS/BAC/dataForSim/',dirName,...
        ';mkdir ',userDir,'/CS/BAC/tmpRuns/',dirName,...
        ';mkdir ',userDir,'/CS/BAC/',dirName])
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%% with correction
  % change this
  correctMeasurementForReadErrorsFlag = 1;
  CorrectReadsNew = 1;
  numProcessors = currNumProcessors;
  basicName = [dirName,'_withCorrection_',num2str(readLength)];
  outFileDirName = [basicName];
  % end change this
  
  basicPrefix = basicName;
  
  %auxData.currDir = [userDir,'/CS/BAC/tmpRuns/',outFileDirName,'/',auxData.tmpFileName];
  %auxData.tmpFileName = outFileDirName;
  %basicSaveName = [auxData.currDir,'/',auxData.tmpFileName,'_dat'];
  
  unix(['mkdir ',userDir,'/CS/tmp/',dirName,'/',outFileDirName,...
        ';mkdir ',userDir,'/CS/BAC/dataForSim/',dirName,'/',outFileDirName,...
        ';mkdir ',userDir,'/CS/BAC/tmpRuns/',dirName,'/',outFileDirName,...
        ';mkdir ',userDir,'/CS/BAC/',dirName,'/',outFileDirName])
  
  auxData = struct;
  
  auxData.readLength = readLength;
  auxData.userDir = userDir;
  auxData.numProcessors = numProcessors;
  
  auxData.correctMeasurementForReadErrorsFlag = correctMeasurementForReadErrorsFlag;
  auxData.CorrectReadsNew = CorrectReadsNew;
  
  auxData.repeatRandomGroups = 10;
  auxData.moveDependentToNextStageFlag = 1;
  
  auxData.tmpFileName = basicName;
  auxData.randomizeSeedAfterCreateReadsFlag = 1;
  auxData.saveName = [userDir,'/CS/BAC/',dirName,'/',outFileDirName,'/',auxData.tmpFileName];
  
  %auxData.currDir = [userDir,'/CS/BAC/tmpRuns/',outFileDirName,'/',auxData.tmpFileName];
  auxData.currDir = [userDir,'/CS/BAC/tmpRuns/',dirName,'/',outFileDirName,'/'];
  
  
  auxData.basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];
  auxData.basicSeqKey= [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
  
  
  
  auxData.brFlag =  1;
  auxData.queueName = 'hour';  
  
  
  auxData.realReadsFromFileFlag = 1;
  
  auxData.reads_data_fileName = [userDir,'/CS/BAC/Peter/',readsDataFile];
    
  
  
  auxData.numBACtoConsider = numBACtoConsider;
  auxData.readLength = readLength;
  auxData.numProcessors = numProcessors;
  auxData.firstFlag = 0;
  auxData.groupSize = 1000;
  auxData.smallestSetOfCollected = 1000;
  auxData.thresholdForCollectingBAC = 10^-3;
  auxDataFileName = [userDir,'/CS/tmp/',dirName,'/',outFileDirName,'/','structFor_',auxData.tmpFileName];
  
  auxData.saveFullDataFlag = 1;
  auxData.upperLimitForRepeat = 60000;
  auxData.repeatWhenLowerThanThisValue = 100000;
  auxData.batchSize = 400;
  
  auxData.inifiniteNumberOfReadsFlag = 0;
  auxData.moveDependentToNextStageFlag = 1;
  
  auxData.createReadsAndQuit=0;
  auxData.loadSavedReadsForLargeReadLengthAndManyReads = 0;
  auxData.dataSetPrimers750Flag = 1;
  %keyboard
  
  save(auxDataFileName,'auxData')
  
  
  auxDataFileName = [userDir,'/CS/tmp/',dirName,'/',outFileDirName,'/','structFor_',basicName];
  currDir = [userDir,'/CS/BAC/tmpRuns/',dirName,'/',outFileDirName]; 
  
  currFileName =        [userDir,'/CS/BAC/dataForSim/',dirName,'/','run_',basicName];
  
  
  curr_fid = fopen(currFileName,'w');
  fprintf(curr_fid,'#!/bin/bash\n');
  fprintf(curr_fid,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
  fprintf(curr_fid,'mkdir $TMPDIR\n');
  fprintf(curr_fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
  
  w1 = [sprintf(['cd ',currDir,';'])];
  fprintf(curr_fid,'%s\n',w1);
  w1 = ['bsub -c 1440 -q week '...
        ' -o ',userDir,'/CS/BAC/dataForSim/',dirName,'/',outFileDirName,'/',basicName,'.o ',...
        ' -e ',userDir,'/CS/BAC/dataForSim/',dirName,'/',outFileDirName,'/',basicName,'.e ',...
        '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/cvx/run_distributeBAC_generalOrFourth.sh  /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b ',auxDataFileName];
  % deleted ' -M 33554432',
  
  fprintf(curr_fid,'%s\n',w1);
  fprintf(curr_fid,'rm -fr $TMPDIR\n');
  fclose(curr_fid);
  clear curr_fid
  
  unix(['chmod 0700 ',currFileName]);
  % end print run file
  %distributeBAC_generalOrFourth(auxDataFileName);

  
  % list file
  writeFile = [userDir,'/CS/BAC/dataForSim/',dirName,'/','listOfR_',basicPrefix];
  fid = fopen(writeFile,'w');
  fprintf(fid,'#!/bin/bash\n');
  fprintf(fid,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
  fprintf(fid,'mkdir $TMPDIR\n');
  fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
  fprintf(fid,'\n');  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%% with corection
  
  fprintf(fid,'%s\n',currFileName);
  
  fprintf(fid,'rm -fr $TMPDIR\n');
  fclose(fid);
  clear fid
  
  
  
  
  
  unix(['chmod 0700 ',writeFile]);
  
end   


