%%%%%%%%%%

clear
auxData = struct;
auxData.outFileDirName = 'canonical200';
auxData.bdf = 0;
auxData.nb = [200];
auxData.pwr = [1]; % power law 1
auxData.readLength = [100];
auxData.numIter = 10;
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
auxData.outFileDirName = 'canonical200_full';
auxData.bdf = 0;
auxData.nb = [200];
auxData.pwr = [1]; % power law 1
auxData.readLength = [100];
auxData.numIter = 100;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10];
auxData.numProcessors{2}([50 100 200 400]) = [50 50 75 150];
auxData.doL1Flag = 1;
runGeneralSplitCorrection(auxData);


%%%
% noise effect




% %%%%%%%%%
% read length

% read length full
clear
auxData = struct;
auxData.outFileDirName = 'VariableReadLength_c200';
auxData.bdf = 0;
auxData.nb = [200];
auxData.pwr = [1]; % power law 1
auxData.readLength = [200 35];
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

clear
auxData = struct;
auxData.outFileDirName = 'VariableReadLength_c200_75_150';
auxData.bdf = 0;
auxData.nb = [200];
auxData.pwr = [1]; % power law 1
auxData.readLength = [75 150];
auxData.numIter = 10;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10];
auxData.numProcessors{2}([35 50 75 150 400]) = [50 50 100 100 150];
auxData.doL1Flag = 1;
runGeneralSplitCorrection(auxData);




clear
auxData = struct;
auxData.outFileDirName = 'VariableReadLength_c200_full_1_20';
auxData.bdf = 0;
auxData.nb = [200];
auxData.pwr = [1]; % power law 1
auxData.readLength = [35 75 150 200];
auxData.numIter = 20;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10];
auxData.numProcessors{2}([35 75 100 150 200]) = [50 100 100 125 150];
auxData.doL1Flag = 1;
runGeneralSplitCorrection(auxData);


delDir('VariableReadLength_c200_full_21_40')

clear
auxData = struct;
auxData.outFileDirName = 'VariableReadLength_c200_full_21_40';
auxData.bdf = 0;
auxData.nb = [200];
auxData.pwr = [1]; % power law 1
auxData.readLength = [35 75 150 200];
auxData.numIter = 20;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10];
auxData.numProcessors{2}([35 75 100 150 200]) = [50 100 100 125 150];
auxData.doL1Flag = 1;

fileNameForFindNumber = 'dataForWithBalance';
auxData.fileNameForFindNumber = fileNameForFindNumber;

auxDataFileNumber = struct;
auxDataFileNumber.maxNumberHourRunning = 450;
auxDataFileNumber.maxNumberHourPending = 500
auxDataFileNumber.maxNumnerWeekRunning = 200; %200
auxDataFileNumber.maxNumberWeekPending = 20; %20
auxDataFileNumber.pauseTimeSeconds = 60;
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
save([userDir,'/CS/BAC/cvx/',fileNameForFindNumber],'auxDataFileNumber');

runGeneralWithLoadBalance(auxData);

%%%%%%%%
clear
auxData = struct;
auxData.outFileDirName = 'VariableReadLength_c200_full_41_100';
auxData.bdf = 0;
auxData.nb = [200];
auxData.pwr = [1]; % power law 1
auxData.readLength = [35 75 150 200];
auxData.numIter = 60;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10];
auxData.numProcessors{2}([35 75 100 150 200]) = [50 100 100 125 150];
auxData.doL1Flag = 1;
auxData.fileNameForFindNumber = 'dataForWithBalance';
runGeneralWithLoadBalance(auxData);


clear
auxData = struct;
auxData.outFileDirName = 'VariableReadLength_c200_50';
auxData.bdf = 0;
auxData.nb = [200];
auxData.pwr = [1]; % power law 1
auxData.readLength = [50];
auxData.numIter = 100;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [50];
auxData.numProcessors{2}([35 50 75 100 150 200]) = [50 50 100 100 125 150];
auxData.doL1Flag = 1;
auxData.fileNameForFindNumber = 'dataForWithBalance';
runGeneralWithLoadBalance(auxData);



% number of bacteria


%%%%%%%%%%%%%%% end number of bacteria


%%%%%%%%%%55
% number of reads 0- should be done again as for less than 10^5 reads - there was no noise added


clear
auxData = struct;
auxData.outFileDirName = 'VariableNumberOfReads_c200';
auxData.bdf = 0;
auxData.nb = [200];
auxData.pwr = [1]; % power law 1
auxData.readLength = [100];
auxData.numIter = 10;
auxData.nr = [10^4 10^5 5*10^5];
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
auxData.outFileDirName = 'VariableNumberOfReads_c200_5x1e4';
auxData.bdf = 0;
auxData.nb = [200];
auxData.pwr = [1]; % power law 1
auxData.readLength = [100];
auxData.numIter = 10;
auxData.nr = [5*10^4];
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
auxData.outFileDirName = 'VariableNumberOfReads_based200_full';
auxData.bdf = 0;
auxData.nb = [200];
auxData.pwr = [1]; % power law 1
auxData.readLength = [100];
auxData.numIter = 100;
auxData.nr = [10^4 2*10^4 5*10^4 10^5 5*10^5];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10];
auxData.numProcessors{2}([50 100 200 400]) = [50 100 75 150];
auxData.doL1Flag = 1;
auxData.fileNameForFindNumber = 'dataForWithBalance';
runGeneralWithLoadBalance(auxData);





% Variable number of Bacetria
clear
auxData = struct;
auxData.outFileDirName = 'VariableNumberOfBacteria_uniform';
auxData.bdf = 0;
auxData.nb = [10 100 200:200:1000];
auxData.pwr = [0]; % uniform
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
auxData.outFileDirName = 'VariableNumberOfBacteria_uniform_full';
auxData.bdf = 0;
auxData.nb = [10 100 200:200:1000];
auxData.pwr = [0]; % uniform
auxData.readLength = [100];
auxData.numIter = 90;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10];
auxData.numProcessors{2}([50 100 200 400]) = [50 50 75 150];
auxData.doL1Flag = 1
auxData.fileNameForFindNumber = 'dataForWithBalance';
runGeneralWithLoadBalance(auxData);


clear
auxData = struct;
auxData.outFileDirName = 'noiseEffect';
auxData.bdf = 0;
auxData.nb = [200];
auxData.pwr = [1]; % power law 1
auxData.readLength = [100];
auxData.numIter = 10;
auxData.nr = [10^5 5*10^5 10^6];
auxData.addNoiseFlagInput = [0 1];
auxData.correctMeasurementForReadErrorsFlagInput = [0 1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10]; %without correction
auxData.numProcessors{2}([50 100 200 400]) = [50 50 75 150];
auxData.doL1Flag = 1;
auxData.fileNameForFindNumber = 'dataForWithBalance';
runGeneralWithLoadBalance(auxData);


clear
auxData = struct;
auxData.outFileDirName = 'noiseEffect2';
auxData.bdf = 0;
auxData.nb = [200];
auxData.pwr = [1]; % power law 1
auxData.readLength = [100];
auxData.numIter = 100;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [0 1];
auxData.correctMeasurementForReadErrorsFlagInput = [0 1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10]; %without correction
auxData.numProcessors{2}([50 100 200 400]) = [50 50 75 150];
auxData.doL1Flag = 1;
auxData.fileNameForFindNumber = 'dataForWithBalance';
runGeneralWithLoadBalance(auxData);

clear
auxData = struct;
auxData.outFileDirName = 'noiseEffect2';
auxData.bdf = 0;
auxData.nb = [200];
auxData.pwr = [1]; % power law 1
auxData.readLength = [100];
auxData.numIter = 100;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [0 1];
auxData.correctMeasurementForReadErrorsFlagInput = [0 1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10]; %without correction
auxData.numProcessors{2}([50 100 200 400]) = [50 50 75 150];
auxData.doL1Flag = 1;
auxData.fileNameForFindNumber = 'dataForWithBalance';
runGeneralWithLoadBalance(auxData);

clear
auxData = struct;
auxData.outFileDirName = 'time_withCorrection';
auxData.bdf = 0;
auxData.nb = [200];
auxData.pwr = [1]; % power law 1
auxData.readLength = [100];
auxData.numIter = 1;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [5]; %without correction
auxData.numProcessors{2}([50 100 200 400]) = [50 50 75 150];
auxData.doL1Flag = 1;
auxData.fileNameForFindNumber = 'dataForWithBalance';
runGeneralWithLoadBalance(auxData);

clear
auxData = struct;
auxData.outFileDirName = 'time_noCorrection';
auxData.bdf = 0;
auxData.nb = [200];
auxData.pwr = [1]; % power law 1
auxData.readLength = [100];
auxData.numIter = 1;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [0];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [5]; %without correction
auxData.numProcessors{2}([50 100 200 400]) = [50 50 75 150];
auxData.doL1Flag = 1;
auxData.fileNameForFindNumber = 'dataForWithBalance';
runGeneralWithLoadBalance(auxData);


clear
auxData = struct;
auxData.outFileDirName = 'similarBACtest';
auxData.bdf = 1;
auxData.nb = [100];
auxData.pwr = [0]; % power law 0
auxData.readLength = [100];
auxData.numIter = 10;
auxData.nr = [10^6];
auxData.addNoiseFlagInput = [1];
auxData.correctMeasurementForReadErrorsFlagInput = [1];
auxData.CorrectReadsNew = 1; % first method
auxData.dataSetPrimers750Flag = 0;
auxData.numProcessors{1} = [10]; %without correction
auxData.numProcessors{2}([50 100 200 400]) = [50 50 75 150];
auxData.doL1Flag = 1;
auxData.fileNameForFindNumber = 'dataForWithBalance';
runGeneralWithLoadBalance(auxData);




%%%%%%%%%%%%%%%%%%%%%
if 1==2
  clear
  % delDir('canonical200Time')
  auxData = struct;
  auxData.outFileDirName = 'canonical200Time';
  auxData.bdf = 0;
  auxData.nb = [200];
  auxData.pwr = [1]; % power law 1
  auxData.readLength = [100];
  auxData.numIter = 1;
  auxData.nr = [10^6];
  auxData.addNoiseFlagInput = [1];
  auxData.correctMeasurementForReadErrorsFlagInput = [1];
  auxData.CorrectReadsNew = 1; % first method
  auxData.dataSetPrimers750Flag = 0;
  auxData.numProcessors{1} = [4];
  auxData.numProcessors{2}([50 100 200 400]) = [4 4 4 4];
  auxData.doL1Flag = 1;

  if ~isfield(auxData,'numProcessors')
    auxData.numProcessors = 10;
  end
  
  defValues = struct;
  defValues.bdf = 0;
  defValues.groupSize = 1000;
  defValues.thresholdForCollectingBAC = 10^-3;
  
  auxData = checkValues(auxData,defValues);

  printStructValues(auxData,'auxData in runGeneral');
  
  userDir = '~/';
  
  unix(['mkdir ',userDir,'/CS/tmp/',auxData.outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/tmpRuns/',auxData.outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/',auxData.outFileDirName])

  prefix = auxData.outFileDirName;
  brFlag = 0;
  
  crGeneralFourth_splitCorrection(prefix,auxData.outFileDirName,auxData.numProcessors,auxData.bdf,auxData.nb,auxData.pwr,auxData.numIter,auxData.nr,auxData.addNoiseFlagInput,auxData.correctMeasurementForReadErrorsFlagInput,auxData.readLength,'hour','week',150000,50000,400,brFlag,1,auxData.thresholdForCollectingBAC,auxData.groupSize,auxData.dataSetPrimers750Flag,auxData.CorrectReadsNew,auxData.doL1Flag);
  crGeneralFourth_splitCorrection(prefix,auxData.outFileDirName,auxData.numProcessors,auxData.bdf,auxData.nb,auxData.pwr,auxData.numIter,auxData.nr,auxData.addNoiseFlagInput,auxData.correctMeasurementForReadErrorsFlagInput,auxData.readLength,'hour','week',150000,50000,400,brFlag,0,auxData.thresholdForCollectingBAC,auxData.groupSize,auxData.dataSetPrimers750Flag,auxData.CorrectReadsNew,auxData.doL1Flag);

  % running - by hand - with correction
  clear
  load ~/CS/tmp/canonical200Time/structFor_canonical200Time_bac_dist_flag_0_Nbac_in_mixture_200_npower_10_readlen_100_numIter_1_Nreads_1000000_Noise_1_Correction_1
  auxData.saveName = '~/CS/BAC/canonical200Time/canonical200Time_bac_dist_flag_0_Nbac_in_mixture_200_npower_10_readlen_100_numIter_1_Nreads_1000000_Noise_1_Correction_1'
  auxData.currDir = '~/CS/BAC/tmpRuns/canonical200Time/canonical200Time_bac_dist_flag_0_Nbac_in_mixture_200_npower_10_readlen_100_numIter_1_Nreads_1000000_Noise_1_Correction_1';
  auxData.reads_data_fileName = '~/CS/BAC/dataForSim/canonical200Time/canonical200Time_bac_dist_flag_0_Nbac_in_mixture_200_npower_10_readlen_100_numIter_1_Nreads_1000000_Noise_1_Correction_1_Nbacmix_200_Nread_1000000_Readlen_100_Npower_1_bacdistflag_0_NoReads';
  auxData.basicSeqNameDir = '~/CS/BAC/datNoNonACGT/packed64/';
  auxData.basicSeqKey = '~/CS/BAC/datNoNonACGT/keyNoNonACGT';
  
  save ~/CS/tmp/canonical200Time/structFor_canonical200Time_bac_dist_flag_0_Nbac_in_mixture_200_npower_10_readlen_100_numIter_1_Nreads_1000000_Noise_1_Correction_1 auxData
  
  t1 = clock;distributeBAC_generalOrFourth(['~/CS/tmp/canonical200Time/structFor_canonical200Time_bac_dist_flag_0_Nbac_in_mixture_200_npower_10_readlen_100_numIter_1_Nreads_1000000_Noise_1_Correction_1']);t2=etime(clock,t1);save ~/time_t1_t2 t1 t2

  
  
end

%%%%%%%%%%%%%%%%%%%

% time - no correction
if 1==2
  clear
  % delDir('canonical200TimeNoCorrection')
  auxData = struct;
  auxData.outFileDirName = 'canonical200TimeNoCorrection_sol';
  auxData.bdf = 0;
  auxData.nb = [200];
  auxData.pwr = [1]; % power law 1
  auxData.readLength = [100];
  auxData.numIter = 1;
  auxData.nr = [10^6];
  auxData.addNoiseFlagInput = [1];
  auxData.correctMeasurementForReadErrorsFlagInput = [0];
  auxData.CorrectReadsNew = 1; % first method
  auxData.dataSetPrimers750Flag = 0;
  %auxData.numProcessors{1} = [4];
  %auxData.numProcessors{2}([50 100 200 400]) = [4 4 4 4];
  auxData.numProcessors{1} = [2];
  auxData.numProcessors{2}([50 100 200 400]) = [2 2 2 2];
  
  auxData.doL1Flag = 1;

  if ~isfield(auxData,'numProcessors')
    auxData.numProcessors = 10;
  end
  
  defValues = struct;
  defValues.bdf = 0;
  defValues.groupSize = 1000;
  defValues.thresholdForCollectingBAC = 10^-3;
  
  auxData = checkValues(auxData,defValues);

  printStructValues(auxData,'auxData in runGeneral');
  
  userDir = '~/';
  
  unix(['mkdir ',userDir,'/CS/tmp/',auxData.outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/tmpRuns/',auxData.outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/',auxData.outFileDirName])

  prefix = auxData.outFileDirName;
  brFlag = 0;
  
  crGeneralFourth_splitCorrection(prefix,auxData.outFileDirName,auxData.numProcessors,auxData.bdf,auxData.nb,auxData.pwr,auxData.numIter,auxData.nr,auxData.addNoiseFlagInput,auxData.correctMeasurementForReadErrorsFlagInput,auxData.readLength,'hour','week',150000,50000,400,brFlag,1,auxData.thresholdForCollectingBAC,auxData.groupSize,auxData.dataSetPrimers750Flag,auxData.CorrectReadsNew,auxData.doL1Flag);
  crGeneralFourth_splitCorrection(prefix,auxData.outFileDirName,auxData.numProcessors,auxData.bdf,auxData.nb,auxData.pwr,auxData.numIter,auxData.nr,auxData.addNoiseFlagInput,auxData.correctMeasurementForReadErrorsFlagInput,auxData.readLength,'hour','week',150000,50000,400,brFlag,0,auxData.thresholdForCollectingBAC,auxData.groupSize,auxData.dataSetPrimers750Flag,auxData.CorrectReadsNew,auxData.doL1Flag);

  % running - by hand - with correction
  clear
  load ~/CS/tmp/canonical200TimeNoCorrection_sol/structFor_canonical200TimeNoCorrection_sol_bac_dist_flag_0_Nbac_in_mixture_200_npower_10_readlen_100_numIter_1_Nreads_1000000_Noise_1_Correction_0
  auxData.saveName = '~/CS/BAC/canonical200TimeNoCorrection_sol/canonical200TimeNoCorrection_sol_bac_dist_flag_0_Nbac_in_mixture_200_npower_10_readlen_100_numIter_1_Nreads_1000000_Noise_1_Correction_0'
  auxData.currDir = '~/CS/BAC/tmpRuns/canonical200TimeNoCorrection_sol/canonical200TimeNoCorrection_sol_bac_dist_flag_0_Nbac_in_mixture_200_npower_10_readlen_100_numIter_1_Nreads_1000000_Noise_1_Correction_0';
  auxData.reads_data_fileName = '~/CS/BAC/dataForSim/canonical200TimeNoCorrection_sol/canonical200TimeNoCorrection_sol_bac_dist_flag_0_Nbac_in_mixture_200_npower_10_readlen_100_numIter_1_Nreads_1000000_Noise_1_Correction_0_Nbacmix_200_Nread_1000000_Readlen_100_Npower_1_bacdistflag_0_NoReads';
  auxData.basicSeqNameDir = '~/CS/BAC/datNoNonACGT/packed64/';
  auxData.basicSeqKey = '~/CS/BAC/datNoNonACGT/keyNoNonACGT';
  
  save ~/CS/tmp/canonical200TimeNoCorrection_sol/structFor_canonical200TimeNoCorrection_sol_bac_dist_flag_0_Nbac_in_mixture_200_npower_10_readlen_100_numIter_1_Nreads_1000000_Noise_1_Correction_0 auxData
  
  t1 = clock;distributeBAC_generalOrFourth(['~/CS/tmp/canonical200TimeNoCorrection_sol/structFor_canonical200TimeNoCorrection_sol_bac_dist_flag_0_Nbac_in_mixture_200_npower_10_readlen_100_numIter_1_Nreads_1000000_Noise_1_Correction_0']);t2=etime(clock,t1);save ~/time_t1_t2_noCorrection_sol t1 t2

  
  
end


if 1==2 %with correction
    clear
  % delDir('canonical200TimeWithCorrection')
  auxData = struct;
  auxData.outFileDirName = 'canonical200TimeWithCorrection_sol';
  auxData.bdf = 0;
  auxData.nb = [200];
  auxData.pwr = [1]; % power law 1
  auxData.readLength = [100];
  auxData.numIter = 1;
  auxData.nr = [10^6];
  auxData.addNoiseFlagInput = [1];
  auxData.correctMeasurementForReadErrorsFlagInput = [1];
  auxData.CorrectReadsNew = 1; % first method
  auxData.dataSetPrimers750Flag = 0;
  %auxData.numProcessors{1} = [4];
  %auxData.numProcessors{2}([50 100 200 400]) = [4 4 4 4];
  auxData.numProcessors{1} = [2];
  auxData.numProcessors{2}([50 100 200 400]) = [2 2 2 2];
  
  auxData.doL1Flag = 1;
  if ~isfield(auxData,'numProcessors')
    auxData.numProcessors = 10;
  end
  
  
  defValues = struct;
  defValues.bdf = 0;
  defValues.groupSize = 1000;
  defValues.thresholdForCollectingBAC = 10^-3;
  
  auxData = checkValues(auxData,defValues);

  printStructValues(auxData,'auxData in runGeneral');
  
  userDir = '~/';
  
  unix(['mkdir ',userDir,'/CS/tmp/',auxData.outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/dataForSim/',auxData.outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/tmpRuns/',auxData.outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/',auxData.outFileDirName])

  prefix = auxData.outFileDirName;
  brFlag = 0;
  
  crGeneralFourth_splitCorrection(prefix,auxData.outFileDirName,auxData.numProcessors,auxData.bdf,auxData.nb,auxData.pwr,auxData.numIter,auxData.nr,auxData.addNoiseFlagInput,auxData.correctMeasurementForReadErrorsFlagInput,auxData.readLength,'hour','week',150000,50000,400,brFlag,1,auxData.thresholdForCollectingBAC,auxData.groupSize,auxData.dataSetPrimers750Flag,auxData.CorrectReadsNew,auxData.doL1Flag);
  crGeneralFourth_splitCorrection(prefix,auxData.outFileDirName,auxData.numProcessors,auxData.bdf,auxData.nb,auxData.pwr,auxData.numIter,auxData.nr,auxData.addNoiseFlagInput,auxData.correctMeasurementForReadErrorsFlagInput,auxData.readLength,'hour','week',150000,50000,400,brFlag,0,auxData.thresholdForCollectingBAC,auxData.groupSize,auxData.dataSetPrimers750Flag,auxData.CorrectReadsNew,auxData.doL1Flag);

  % running - by hand - with correction
  clear
  load ~/CS/tmp/canonical200TimeWithCorrection_sol/structFor_canonical200TimeWithCorrection_sol_bac_dist_flag_0_Nbac_in_mixture_200_npower_10_readlen_100_numIter_1_Nreads_1000000_Noise_1_Correction_1
  auxData.saveName = '~/CS/BAC/canonical200TimeWithCorrection_sol/canonical200TimeWithCorrection_sol_bac_dist_flag_0_Nbac_in_mixture_200_npower_10_readlen_100_numIter_1_Nreads_1000000_Noise_1_Correction_1'
  auxData.currDir = '~/CS/BAC/tmpRuns/canonical200TimeWithCorrection_sol/canonical200TimeWithCorrection_sol_bac_dist_flag_0_Nbac_in_mixture_200_npower_10_readlen_100_numIter_1_Nreads_1000000_Noise_1_Correction_1';
  auxData.reads_data_fileName = '~/CS/BAC/dataForSim/canonical200TimeWithCorrection_sol/canonical200TimeWithCorrection_sol_bac_dist_flag_0_Nbac_in_mixture_200_npower_10_readlen_100_numIter_1_Nreads_1000000_Noise_1_Correction_1_Nbacmix_200_Nread_1000000_Readlen_100_Npower_1_bacdistflag_0_NoReads';
  auxData.basicSeqNameDir = '~/CS/BAC/datNoNonACGT/packed64/';
  auxData.basicSeqKey = '~/CS/BAC/datNoNonACGT/keyNoNonACGT';
  
  save ~/CS/tmp/canonical200TimeWithCorrection_sol/structFor_canonical200TimeWithCorrection_sol_bac_dist_flag_0_Nbac_in_mixture_200_npower_10_readlen_100_numIter_1_Nreads_1000000_Noise_1_Correction_1 auxData
  
  t1 = clock;distributeBAC_generalOrFourth(['~/CS/tmp/canonical200TimeWithCorrection_sol/structFor_canonical200TimeWithCorrection_sol_bac_dist_flag_0_Nbac_in_mixture_200_npower_10_readlen_100_numIter_1_Nreads_1000000_Noise_1_Correction_1']);t2=etime(clock,t1);save ~/time_t1_t2_withCorrection_sol t1 t2

  
end