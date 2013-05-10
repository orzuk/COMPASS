%%%%%%%

%keyboard
if exist('found')
  %unix(['rm -f ',auxData.currDir,'/',auxData.tmpFileName,'_*']);
  unix(['rm -fr ',auxData.currDir,'/']);
  %unix(['rm -f ',auxData.currDir,'/run_',auxData.tmpFileName,'_*']);
end


if 1==2 % large
  clear
  userDir = getuserdir;
  outFileDirName = 'test';
  basicName = outFileDirName;
  outfilename = [userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',basicName];
  unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName,...
       ';mkdir ',userDir,'/CS/BAC/',outFileDirName])
  
  load([userDir,'/CS/BAC/cellsimbac2'])
  Nbac_in_mixture = 100;
  Nreads = Inf;
  readlen = 50;
  npower = 0;
  bac_dist_flag = 0;
  
  N = 410849;
  createReadsFlag = 0;
  [strinfilename] = fun_prep_data_for_simulations_loadSpecificBAC2([],[],indsimbac_vec,Nbac_in_mixture,Nreads,readlen,npower,bac_dist_flag,[],outfilename,N,createReadsFlag);
  
  
  
  
  auxData = struct;

  auxData.repeatRandomGroups = 10;
  auxData.moveDependentToNextStageFlag = 1;
  
  auxData.tmpFileName = outFileDirName;
  auxData.randomizeSeedAfterCreateReadsFlag = 1;
  auxData.saveName = [userDir,'/CS/BAC/',outFileDirName,'/',auxData.tmpFileName];
  auxData.inifiniteNumberOfReadsFlag = 1;
  auxData.currDir = [userDir,'/CS/BAC/tmpRuns/',outFileDirName,'/',auxData.tmpFileName];
  auxData.basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];
  auxData.basicSeqKey= [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
  auxData.createReadsFlag = 1;
  auxData.brFlag =  0;
  auxData.queueName = '';  
  if auxData.brFlag
    auxData.reads_data_fileName = strinfilename;
  else
    tmpStr = strrep(strinfilename,userDir,userDir);
    auxData.reads_data_fileName = tmpStr;
  end
  
  auxData.readLength = readlen;
  auxData.numProcessors = 6;
  auxData.firstFlag = 0;
  auxData.groupSize = 1000;
  auxData.smallestSetOfCollected = 1000;
  auxData.numBACtoConsider = N;
  auxData.thresholdForCollectingBAC = 1e-3;
  auxDataFileName = [userDir,'/CS/tmp/',outFileDirName,'/','structFor_',auxData.tmpFileName];
  auxData.addNoiseFlag = 0;
  auxData.correctMeasurementForReadErrorsFlag = 0;
  auxData.saveFullDataFlag = 1;
  auxData.upperLimitForRepeat = 150000;
  auxData.repeatWhenLowerThanThisValue = 50000;
  auxData.batchSize = 400;
  
  save(auxDataFileName,'auxData')
  
  distributeBAC_generalOrFourth(auxDataFileName);
  
end

if 1==2 %small
  clear
  userDir = getuserdir;
  outFileDirName = 'ffSmall';
  basicName = outFileDirName;
  outfilename = [userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',basicName];
  unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,...
        ';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,...
        ';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName])
  
  load([userDir,'/CS/BAC/cellsimbac2'])
  Nbac_in_mixture = 10;
  Nreads = Inf;
  
  
  readlen = 100;
  npower = 0;
  bac_dist_flag = 0;
  
  %N = 10000;
  N = 2000;
  createReadsFlag = 0;
  [strinfilename] = fun_prep_data_for_simulations_loadSpecificBAC2([],[],indsimbac_vec,Nbac_in_mixture,Nreads,readlen,npower,bac_dist_flag,[],outfilename,N,createReadsFlag);
  
  
  
  
  auxData = struct;

  auxData.repeatRandomGroups = 10;
  auxData.moveDependentToNextStageFlag = 1;
  
  auxData.tmpFileName = outFileDirName;
  auxData.randomizeSeedAfterCreateReadsFlag = 1;
  auxData.saveName = [userDir,'/CS/BAC/',outFileDirName,'/',auxData.tmpFileName];
  auxData.inifiniteNumberOfReadsFlag = 1;
  auxData.currDir = [userDir,'/CS/BAC/tmpRuns/',outFileDirName,'/',auxData.tmpFileName];
  auxData.basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];
  auxData.basicSeqKey= [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
  auxData.createReadsFlag = 1;
  auxData.brFlag =  0;
  auxData.queueName = '';  
  if auxData.brFlag
    auxData.reads_data_fileName = strinfilename;
  else
    tmpStr = strrep(strinfilename,userDir,userDir);
    auxData.reads_data_fileName = tmpStr;
  end

  auxData.readLength = 100;
  auxData.numProcessors = 2;
  auxData.firstFlag = 0;
  auxData.groupSize = 1000;
  auxData.smallestSetOfCollected = 1000;
  auxData.numBACtoConsider = N;
  auxData.thresholdForCollectingBAC = 1e-3;
  auxDataFileName = [userDir,'/CS/tmp/',outFileDirName,'/','structFor_',auxData.tmpFileName];
  
  auxData.addNoiseFlag = 0;
  auxData.correctMeasurementForReadErrorsFlag = 0;
  auxData.saveFullDataFlag = 1;
  auxData.upperLimitForRepeat = 6000;
  auxData.repeatWhenLowerThanThisValue = 2000;
  auxData.batchSize = 20;
  
  auxData.doL1Flag = 1;
  auxData.dataSetPrimers750Flag = 0;
  save(auxDataFileName,'auxData')
  
  distributeBAC_generalOrFourth(auxDataFileName);
  
end

if 1==2 % large with 1000 bacteria group of 5000
  clear
  profile -memory on
  userDir = getuserdir;
  outFileDirName = 'test1000';
  basicName = outFileDirName;
  outfilename = [userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',basicName];
  unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName,...
       ';mkdir ',userDir,'/CS/BAC/',outFileDirName])
  
  load([userDir,'/CS/BAC/cellsimbac2'])
  Nbac_in_mixture = 1000;
  Nreads = Inf;
  readlen = 400;
  npower = 0;
  bac_dist_flag = 0;
  
  N = 410849;
  createReadsFlag = 0;
  [strinfilename] = fun_prep_data_for_simulations_loadSpecificBAC2([],[],indsimbac_vec,Nbac_in_mixture,Nreads,readlen,npower,bac_dist_flag,[],outfilename,N,createReadsFlag);
  
  
  
  
  auxData = struct;

  auxData.repeatRandomGroups = 10;
  auxData.moveDependentToNextStageFlag = 1;
  
  auxData.tmpFileName = outFileDirName;
  auxData.randomizeSeedAfterCreateReadsFlag = 1;
  auxData.saveName = [userDir,'/CS/BAC/',outFileDirName,'/',auxData.tmpFileName];
  auxData.inifiniteNumberOfReadsFlag = 1;
  auxData.currDir = [userDir,'/CS/BAC/tmpRuns/',outFileDirName,'/',auxData.tmpFileName];
  auxData.basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];
  auxData.basicSeqKey= [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
  auxData.createReadsFlag = 1;
  auxData.brFlag =  0;
  auxData.queueName = '';  
  if auxData.brFlag
    auxData.reads_data_fileName = strinfilename;
  else
    tmpStr = strrep(strinfilename,userDir,userDir);
    auxData.reads_data_fileName = tmpStr;
  end
  
  auxData.readLength = readlen;
  auxData.numProcessors = 6;
  auxData.firstFlag = 0;
  auxData.groupSize = 5000;
  auxData.smallestSetOfCollected = 1000;
  auxData.numBACtoConsider = N;
  auxData.thresholdForCollectingBAC = 1e-4;
  auxDataFileName = [userDir,'/CS/tmp/',outFileDirName,'/','structFor_',auxData.tmpFileName];
  auxData.addNoiseFlag = 0;
  auxData.correctMeasurementForReadErrorsFlag = 0;
  auxData.saveFullDataFlag = 1;
  auxData.upperLimitForRepeat = 150000;
  auxData.repeatWhenLowerThanThisValue = 50000;
  auxData.batchSize = 400/5;
  
  save(auxDataFileName,'auxData')
  
  distributeBAC_generalOrFourth(auxDataFileName);
  
end


if 1==2 % finite reads long read length
  clear
  
  userDir = getuserdir;
  outFileDirName = 'testFiniteReads';
  basicName = outFileDirName;
  outfilename = [userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',basicName];
  unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName,...
       ';mkdir ',userDir,'/CS/BAC/',outFileDirName])
  
  load([userDir,'/CS/BAC/cellsimbac2'])
  Nbac_in_mixture = 100;
  Nreads = 10^5;
  readlen = 400;
  npower = 0;
  bac_dist_flag = 0;
  
  N = 410849;
  createReadsFlag = 0;
  [strinfilename] = fun_prep_data_for_simulations_loadSpecificBAC2([],[],indsimbac_vec,Nbac_in_mixture,Nreads,readlen,npower,bac_dist_flag,[],outfilename,N,createReadsFlag);
  
  
  
  
  auxData = struct;

  auxData.repeatRandomGroups = 10;
  auxData.moveDependentToNextStageFlag = 1;
  
  auxData.tmpFileName = outFileDirName;
  auxData.randomizeSeedAfterCreateReadsFlag = 1;
  auxData.saveName = [userDir,'/CS/BAC/',outFileDirName,'/',auxData.tmpFileName];
  
  auxData.currDir = [userDir,'/CS/BAC/tmpRuns/',outFileDirName,'/',auxData.tmpFileName];
  auxData.basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];
  auxData.basicSeqKey= [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
  auxData.createReadsFlag = 1;
  auxData.brFlag =  0;
  auxData.queueName = '';  
  if auxData.brFlag
    auxData.reads_data_fileName = strinfilename;
  else
    tmpStr = strrep(strinfilename,userDir,userDir);
    auxData.reads_data_fileName = tmpStr;
  end
  
  auxData.readLength = readlen;
  auxData.numProcessors = 6;
  auxData.firstFlag = 0;
  auxData.groupSize = 1000;
  auxData.smallestSetOfCollected = 1000;
  auxData.numBACtoConsider = N;
  auxData.thresholdForCollectingBAC = 10^-3;
  auxDataFileName = [userDir,'/CS/tmp/',outFileDirName,'/','structFor_',auxData.tmpFileName];
  
  auxData.saveFullDataFlag = 1;
  auxData.upperLimitForRepeat = 150000;
  auxData.repeatWhenLowerThanThisValue = 50000;
  auxData.batchSize = 400;
  
  
  ErrorStruct = []; 
  ErrorStruct.error_model = 'exponential';
  ErrorStruct.baseline_error = 0.005; 
  ErrorStruct.final_error = 0.03; 
  
  p_one_nuc_error_mat = [ -1   0.3  0.22  0.18
                      0.5    -1  0.22   0.6
                      0.35  0.15    -1  0.22
                      0.15  0.55  0.56   -1]; % , 16, length(p));
  
  ErrorStruct.p_one_nuc_error_mat = p_one_nuc_error_mat;
  auxData.ErrorStruct = ErrorStruct;
  auxData.addNoiseFlag = 1;
  auxData.correctMeasurementForReadErrorsFlag = 0;
  auxData.inifiniteNumberOfReadsFlag = 0;
  
  save(auxDataFileName,'auxData')
  
  distributeBAC_generalOrFourth(auxDataFileName);
  
end



if 1==2 % 10^6 reads long read length
  clear
  %profile -memory on
  userDir = getuserdir;
  outFileDirName = 'testNoCorrection';
  basicName = outFileDirName;
  outfilename = [userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',basicName];
  unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName,...
       ';mkdir ',userDir,'/CS/BAC/',outFileDirName])
  
  load([userDir,'/CS/BAC/cellsimbac2'])
  Nbac_in_mixture = 1000;
  Nreads = 10^6;
  readlen = 100;
  npower = 0;
  bac_dist_flag = 0;
  
  N = 3000%410849;
  createReadsFlag = 0;
  [strinfilename] = fun_prep_data_for_simulations_loadSpecificBAC2([],[],indsimbac_vec,Nbac_in_mixture,Nreads,readlen,npower,bac_dist_flag,[],outfilename,N,createReadsFlag);
  
  
  
  
  auxData = struct;
  auxData.userDir = userDir;
  
  auxData.repeatRandomGroups = 2;
  auxData.moveDependentToNextStageFlag = 1;
  
  auxData.tmpFileName = outFileDirName;
  auxData.randomizeSeedAfterCreateReadsFlag = 1;
  auxData.saveName = [userDir,'/CS/BAC/',outFileDirName,'/',auxData.tmpFileName];
  
  auxData.currDir = [userDir,'/CS/BAC/tmpRuns/',outFileDirName,'/',auxData.tmpFileName];
  auxData.basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];
  auxData.basicSeqKey= [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
  auxData.createReadsFlag = 1;
  auxData.brFlag =  0;
  auxData.queueName = '';  
  if auxData.brFlag
    auxData.reads_data_fileName = strinfilename;
  else
    tmpStr = strrep(strinfilename,userDir,userDir);
    auxData.reads_data_fileName = tmpStr;
  end
  
  auxData.readLength = readlen;
  auxData.numProcessors = 1;
  auxData.firstFlag = 0;
  auxData.groupSize = 1000;
  auxData.smallestSetOfCollected = 1000;
  auxData.numBACtoConsider = N;
  auxData.thresholdForCollectingBAC = 10^-3;
  auxDataFileName = [userDir,'/CS/tmp/',outFileDirName,'/','structFor_',auxData.tmpFileName];
  
  auxData.saveFullDataFlag = 1;
  auxData.upperLimitForRepeat = 60000;
  auxData.repeatWhenLowerThanThisValue = 50000;
  auxData.batchSize = 400;
  
  
  ErrorStruct = []; 
  ErrorStruct.error_model = 'exponential';
  ErrorStruct.baseline_error = 0.005; 
  ErrorStruct.final_error = 0.03; 
  
  p_one_nuc_error_mat = [ -1   0.3  0.22  0.18
                      0.5    -1  0.22   0.6
                      0.35  0.15    -1  0.22
                      0.15  0.55  0.56   -1]; % , 16, length(p));
  
  ErrorStruct.p_one_nuc_error_mat = p_one_nuc_error_mat;
  auxData.ErrorStruct = ErrorStruct;
  auxData.addNoiseFlag = 1;
  auxData.inifiniteNumberOfReadsFlag = 0;

  %%% new
  auxData.correctMeasurementForReadErrorsFlag = 1;
  auxData.CorrectReadsNew = 1;
  % end new
  
  auxData.createReadsAndQuit=0;
  auxData.loadSavedReadsForLargeReadLengthAndManyReads = 0;
  auxData.dataSetPrimers750Flag = 0;
  
  save(auxDataFileName,'auxData')
  
  distributeBAC_generalOrFourth(auxDataFileName);
  
end

%%%%%%%%%%%%%%555
if 1==2 %smallBr
  clear
  userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
  outFileDirName = 'ffSmall';
  basicName = outFileDirName;
  outfilename = [userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',basicName];
  unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,...
        ';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,...
        ';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName])
  
  load([userDir,'/CS/BAC/cellsimbac2'])
  Nbac_in_mixture = 10;
  Nreads = Inf;
  
  
  readlen = 100;
  npower = 0;
  bac_dist_flag = 0;
  
  %N = 10000;
  N = 2000;
  createReadsFlag = 0;
  [strinfilename] = fun_prep_data_for_simulations_loadSpecificBAC2([],[],indsimbac_vec,Nbac_in_mixture,Nreads,readlen,npower,bac_dist_flag,[],outfilename,N,createReadsFlag);
  
  
  
  
  auxData = struct;

  auxData.repeatRandomGroups = 10;
  auxData.moveDependentToNextStageFlag = 1;
  
  auxData.tmpFileName = outFileDirName;
  auxData.randomizeSeedAfterCreateReadsFlag = 1;
  auxData.saveName = [userDir,'/CS/BAC/',outFileDirName,'/',auxData.tmpFileName];
  auxData.inifiniteNumberOfReadsFlag = 1;
  auxData.currDir = [userDir,'/CS/BAC/tmpRuns/',outFileDirName,'/',auxData.tmpFileName];
  auxData.basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];
  auxData.basicSeqKey= [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
  auxData.createReadsFlag = 1;
  auxData.brFlag =  1;
  auxData.queueName = '';  
  if auxData.brFlag
    auxData.reads_data_fileName = strinfilename;
  else
    tmpStr = strrep(strinfilename,userDir,userDir);
    auxData.reads_data_fileName = tmpStr;
  end

  auxData.readLength = 100;
  auxData.numProcessors = 1;
  auxData.firstFlag = 0;
  auxData.groupSize = 1000;
  auxData.smallestSetOfCollected = 1000;
  auxData.numBACtoConsider = N;
  auxData.thresholdForCollectingBAC = 1e-3;
  auxDataFileName = [userDir,'/CS/tmp/',outFileDirName,'/','structFor_',auxData.tmpFileName];
  
  auxData.addNoiseFlag = 0;
  auxData.correctMeasurementForReadErrorsFlag = 0;
  auxData.saveFullDataFlag = 1;
  auxData.upperLimitForRepeat = 6000;
  auxData.repeatWhenLowerThanThisValue = 2000;
  auxData.batchSize = 20;
  
  auxData.doL1Flag = 1;
  auxData.dataSetPrimers750Flag = 0;
  save(auxDataFileName,'auxData')
  
  distributeBAC_generalOrFourth(auxDataFileName);
  
end
