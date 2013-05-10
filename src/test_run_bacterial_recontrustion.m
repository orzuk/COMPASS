% Simulate bacterias, run reconstruction

% In the broad: userDir = /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/'


% Run a small example
clear
AssignGeneralConstants;
switch machine
    case UNIX
        userDir = getuserdir;
        userDir = fullfile(userDir, 'CS');
    case PC
        userDir = '../../compressed_sensing/metagenomics/nextgen/data/';
        
end
profile -memory on

outFileDirName = 'testNoCorrection';
basicName = outFileDirName;
outfilename = [userDir,'/BAC/dataForSim/',outFileDirName,'/',basicName];
unix(['mkdir ',userDir,'/tmp/',outFileDirName,...
    ';mkdir ',userDir,'/BAC/dataForSim/',outFileDirName,...
    ';mkdir ',userDir,'/BAC/tmpRuns/',outFileDirName,...
    ';mkdir ',userDir,'/BAC/',outFileDirName])
load([userDir,'/BAC/cellsimbac2'])
Nbac_in_mixture = 1000;
Nreads = 10^6;
readlen = 50;
npower = 0;
bac_dist_flag = 0;

N = 3000%410849;
createReadsFlag = 0; % don't make reads
[strinfilename] = fun_prep_data_for_simulations_loadSpecificBAC2( ...
    [],[],indsimbac_vec,Nbac_in_mixture,Nreads,readlen,npower,bac_dist_flag, ...
    [],outfilename,N,createReadsFlag); % simulate reads



%%%%%%%%%%%%%%%  Set Algorithm Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
auxData = struct;

algorithm_str_vec = {'Noam', 'L2', 'Greedy'} % loop on different algorithms

for i = 1:length(algorithm_str_vec) % loop on algorithm
    auxData.algorithm = algorithm_str_vec{i};
    
    switch auxData.algorithm
        case 'Noam' % currently used algorithm
            auxData.userDir = userDir;
            
            auxData.repeatRandomGroups = 2;
            auxData.moveDependentToNextStageFlag = 1;
            
            auxData.tmpFileName = outFileDirName;
            auxData.randomizeSeedAfterCreateReadsFlag = 1;
            auxData.saveName = fullfile(userDir,'/BAC/', outFileDirName, auxData.tmpFileName);
            
            auxData.currDir = fullfile(userDir, '/BAC/tmpRuns/', outFileDirName, auxData.tmpFileName);
            auxData.basicSeqNameDir = fullfile(userDir,'/BAC/datNoNonACGT/packed64/');
            auxData.basicSeqKey= fullfile(userDir,'/BAC/datNoNonACGT/keyNoNonACGT');
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
            auxDataFileName = fullfile(userDir,'CS/tmp',outFileDirName, ['structFor_',auxData.tmpFileName]);
            
            auxData.saveFullDataFlag = 1;
            auxData.upperLimitForRepeat = 60000;
            auxData.repeatWhenLowerThanThisValue = 50000;
            auxData.batchSize = 400;
            
        case 'Greedy' % use greedy algorithm
            
            
        case 'L2' % use L2 algorithm
    end % swith algorithm
    
    %%%%%%%%%%%%%%%%%%%%%% Set Read Error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ReadErrorStruct = [];
    ReadErrorStruct.error_model = 'exponential';
    ReadErrorStruct.baseline_error = 0.005;
    ReadErrorStruct.final_error = 0.03;
    p_one_nuc_error_mat = [ -1   0.3  0.22  0.18
        0.5    -1  0.22   0.6
        0.35  0.15    -1  0.22
        0.15  0.55  0.56   -1]; % Errors chosen from ???
    ReadErrorStruct.p_one_nuc_error_mat = p_one_nuc_error_mat;
    auxData.ReadErrorStruct = ReadErrorStruct;
    auxData.addNoiseFlag = 1;
    auxData.inifiniteNumberOfReadsFlag = 0;
    %%%%%%%%%%%%%%%%%%%%%% End Set Read Error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%% new
    auxData.correctMeasurementForReadErrorsFlag = 1;
    auxData.CorrectReadsNew = 1;
    % end new
    
    auxData.createReadsAndQuit=0;
    auxData.loadSavedReadsForLargeReadLengthAndManyReads = 0;
    auxData.dataSetPrimers750Flag = 0;
    
    my_mkdir(dir_from_file_name(auxDataFileName));
    save(auxDataFileName,'auxData')
    
    distributeBAC_generalOrFourth(auxDataFileName); % Run Reconstruction (distribute work to several cpu's). Does this part also simulate reads??
end % loop on different algorithms

% Evaluate Results: Compute Reconstruction Error
% Here need code from Amnon !!!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
