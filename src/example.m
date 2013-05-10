addpath(genpath([userDir,'/CS/BAC/cvx/']))


%%%%%%%%%%
% simulation of real reads in fasta file
clear

numBacteriaInDataBase = 410849;
bacteria_in_mix = [1 10 1000];
correctWeight(bacteria_in_mix) = 1/length(bacteria_in_mix);

% read errrors - cite paper
ErrorStruct = []; ErrorStruct.error_model = 'exponential';
ErrorStruct.baseline_error = 0.005; 
ErrorStruct.final_error = 0.03; 

p_one_nuc_error_mat = [ -1   0.3  0.22  0.18
                    0.5    -1  0.22   0.6
                    0.35  0.15    -1  0.22
                    0.15  0.55  0.56   -1]; 

ErrorStruct.p_one_nuc_error_mat = p_one_nuc_error_mat;


% create reads


basicDir = '~/CS/mFilesBAC/package/'; % change to specific location
auxData = struct;
auxData.correctWeight = correctWeight; % vector 
auxData.Nreads = 10^5;
auxData.readLength = 100;
auxData.seed_in = 1;% seed for reads creation
auxData.basicSeqNameDir = [basicDir,'database/datNoNonACGT/packed64/']
auxData.basicSeqKey= [basicDir,'database/datNoNonACGT/keyNoNonACGT']
auxData.addNoiseFlag = 1; % 1 for adding noise, 0= no noise
auxData.ErrorStruct = ErrorStruct;

% create reads in a packed format
[uniqueReads,uniqueReads_length] = createReads_package(auxData);

% for real - start from this point on

% pack the reads and then

% SHRIMP using one core - without correction
% or parfor if parallel is available


% parameters:

setParameters = struct;
setParameters.numBacteriaInDataBase = numBacteriaInDataBase;
setParameters.thresholdForCollectingBAC = 1e-3;
setParameters.upperLimitForRepeat = 150000;
setParameters.repeatWhenLowerThanThisValue = 20000;
setParameters.smallestSetOfCollected = 1000;% after collecting the bacteria from each group - continue iterating id their number is smaller than smallestSetOfCollected


