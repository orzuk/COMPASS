function distributeBAC_general(auxDataFileName)
% reads are red
% can load specific reads - just like needed in similar BAC! just make sure the reads are called red

% create stand alone run files and then run them
%cd ~/CS; addpath('~/CS/mFilesBAC'); mcc -m distributeBAC_general  -d /homes/csfaculty/shental/CS/BAC/cvx  -a currReads.m -a findParametersReads.m -a fun_prep_data_for_simulations_loadSpecificBAC2.m -a get_duplicates.m -a get_duplicates2.m -a getuserdir.m -a iterateParallelDistributedSequenceFiles2.m  -a myHistForBAC.m -a myUniqueUINT8.m -a randi.m -a newReadsBasedSeed2.m -a prepareGroupOf1000DistributedSequenceFiles2.m -a runOneGroupOf1000post.m -a runOneNodeCompiled3.m -a runOneGroupOf1000ForCompilation.m

%cd ~/CS; addpath('~/CS/mFilesBAC'); mcc -m distributeBAC_general  -d /ph2users/fenoam/CS/BAC/cvx  -a currReads.m -a findParametersReads.m -a fun_prep_data_for_simulations_loadSpecificBAC2.m -a get_duplicates.m -a get_duplicates2.m -a getuserdir.m -a iterateParallelDistributedSequenceFiles2.m  -a myHistForBAC.m -a myUniqueUINT8.m -a randi.m -a newReadsBasedSeed2.m -a prepareGroupOf1000DistributedSequenceFiles2.m -a runOneGroupOf1000post.m -a runOneNodeCompiled3.m -a runOneGroupOf1000ForCompilation.m

%keyboard

if isdeployed  %load a file which has a variable by the name auxData 
  load(auxDataFileName)
else
  load(auxDataFileName)
end


if ~isfield(auxData,'keepOriginalOrderFlag') % keepOriginalOrderFlag for  checking in the original order
  auxData.keepOriginalOrderFlag = 0;
end


if ~isfield(auxData,'randomizeSeedAfterCreateReadsFlag') % randomizeSeedAfterCreateReadsFlag for considering a certain bacteria from a group
  auxData.randomizeSeedAfterCreateReadsFlag = 0;
end

% add to auxData stuff
if ~isfield(auxData,'groupSize') % partition into groups of groupSize bacteria
  auxData.groupSize = 1000;
end

if ~isfield(auxData,'smallestSetOfCollected') % after collecting the bacteria from each group - continue iterating id their number is smaller than smallestSetOfCollected
  auxData.smallestSetOfCollected = 1000;
end

if ~isfield(auxData,'thresholdForCollectingBAC') % thresholdForCollectingBAC for considering a certain bacteria from a group
  auxData.thresholdForCollectingBAC = 1e-3;
end

if ~isfield(auxData,'numBACtoConsider')
  auxData.numBACtoConsider = 455055;
end

if ~isfield(auxData,'createReadsFlag')
  auxData.createReadsFlag = 0;
end

%%%%%% %end add to auxData stuff


if ~isfield(auxData,'brFlag')
  auxData.brFlag = 0;
end

if ~isfield(auxData,'queueName')
  auxData.queueName = [];
end

if ~isfield(auxData,'firstFlag')
  auxData.firstFlag = 0;
end


if auxData.brFlag
  userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
  disp('run in Br')
  parallelType = 'cluster';
else
  userDir = getuserdir;
  if ~isempty(findstr(userDir,'fenoam')) %WIS
    parallelType = 'cluster';
    disp('runs in WIS')
  else %OU
    parallelType = 'local';
    disp('run in OU')
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end load variables

fn = fieldnames(auxData);
fprintf('%%%%%%%%%%%\n');
fprintf('data:\n');
for i=1:length(fn)
  if isnumeric(auxData.(fn{i}))
    fprintf('%s: %d\n',fn{i},auxData.(fn{i}));
  else
    fprintf('%s: %s\n',fn{i},auxData.(fn{i}));
  end
end
fprintf('end data\n')
fprintf('%%%%%%%%%%%%\n')

%%%%%%%%%%%%%%%%%%555

if ~isdeployed && auxData.firstFlag==1 
  addpath([userDir,'/CS/mFilesBAC/']);
  addpath(genpath([userDir,'/CS/BAC/cvx/']))
end


%keyboard

%unix(['mkdir ',auxData.currDir]);
runFileName = [auxData.currDir,'/run_',auxData.tmpFileName];
basicSaveName = [auxData.currDir,'/',auxData.tmpFileName,'_dat'];
basicSeqNameDir = [userDir,'/CS/BAC/dat450000/'];
basicSeqKey =  [userDir,'/CS/BAC/dat450000/key450000'];

% create the reads
if auxData.createReadsFlag 
  [junk,red] = newReadsBasedSeed2(auxData.reads_data_fileName,userDir,basicSeqNameDir,basicSeqKey,auxData.numBACtoConsider);
else % load the reads themselves
  load(auxData.reads_data_fileName)
end
disp('assumes that red is uint8');

if ~strcmp(class(red),'uint8')  % Check if data is a uint8
  red = uint8(red);         
  disp('converted red to uint8- should have been saved are uint8!!! in distributeBAC_general.m')
end

%keyboard
[uniqueReads,uniqueReads_inds] = myUniqueUINT8(red);
clear red

uniqueReads_length = zeros(size(uniqueReads,1),1);
for i=1:size(uniqueReads_length,1)
  uniqueReads_length(i) = length(uniqueReads_inds{i});
end
clear uniqueReads_inds

% for distributed - using distributed sequence files
%keyboard

if auxData.randomizeSeedAfterCreateReadsFlag
  disp('randomized state, after the reads were created')
  rand('seed',sum(100*clock));
else
  disp('use the original seed')
end

%keyboard
clear store* 
k = 1;
curr_kp = 1:auxData.numBACtoConsider;
while length(curr_kp)>auxData.smallestSetOfCollected
  [currX,currSumRelevantReads] = iterateParallelDistributedSequenceFiles2(uniqueReads,uniqueReads_length,curr_kp,basicSeqNameDir,basicSeqKey,parallelType,[basicSaveName,'_k_',num2str(k)],userDir,[runFileName,'_k_',num2str(k)],auxData);
  next_kp = curr_kp(find(currX>auxData.thresholdForCollectingBAC));
  store_kp{k} = curr_kp;
  store_X{k} = currX;
  store_sumRelevantReads{k} = currSumRelevantReads;
  k = k +1;
  
  curr_kp = next_kp;
end
store_kp{k} = curr_kp;
store_X{k} = [];

% solve again
lengths_kp = zeros(1,length(store_kp));
for i=1:length(lengths_kp)
  lengths_kp(i) = length(store_kp{i});
end
%1
%keyboard
do = find(lengths_kp<auxData.smallestSetOfCollected*2);
saveSolveAgainDataName = [auxData.currDir,'/',auxData.tmpFileName,'_solveAgain'];
save(saveSolveAgainDataName,'store_kp','store_X','store_sumRelevantReads','basicSeqNameDir','basicSeqKey','uniqueReads','uniqueReads_length','do','auxData');

%pause
% deleted since contained in auxData: 'n','tmpFileName','createReadsFlag'

if strcmp(parallelType,'local')
  solveAgain2(saveSolveAgainDataName,auxData.saveName);
else
  if ~isempty(findstr(userDir,'/storage/')) % WICC
    unix([userDir,'/CS/BAC/cvx/run_solveAgain2_WIS.sh /storage/MATLAB/R2010b/ ',saveSolveAgainDataName,' ',auxData.saveName]);
    
  elseif ~isempty(findstr(userDir,'/seq/'))  % Br
    unix([userDir,'/CS/BAC/cvx/run_solveAgain2_Br.sh /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b ',saveSolveAgainDataName,' ',auxData.saveName]);
  elseif ~isempty(findstr(userDir,'/shental/'))
    unix([userDir,'/CS/BAC/cvx/run_solveAgain2_OU.sh /u/01/LARGE_PKG/MATLAB_2010B ',saveSolveAgainDataName,' ',auxData.saveName]);
  end

end

%quit




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 1==2
  
  clear
  load([getuserdir,'/CS/BAC/cellsimbac2'])
  Nbac_in_mixture = 100;
  Nreads = 2e6;
  readlen = 50;
  npower = 0.5;
  bac_dist_flag = 1;
  outfilename = '~/CS/BAC/dataForSim/simSec';
  N = 455055;
  createReadsFlag = 0;
  [strinfilename] = fun_prep_data_for_simulations_loadSpecificBAC2([],[],indsimbac_vec,Nbac_in_mixture,Nreads,readlen,npower,bac_dist_flag,[],outfilename,N,createReadsFlag);
  
  
  % OU
  clear
  addpath('~/CS/mFilesBAC')
  
  auxData = struct;
  auxData.reads_data_fileName = '~/CS/BAC/dataForSim/simSec_Nbacmix_100_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads';
  auxData.currDir = '/homes/csfaculty/shental/CS/BAC/tmpRuns/testNew';
  auxData.readLength = 50;
  auxData.tmpFileName = 'testNew';
  auxData.saveName = '~/CS/BAC/testNew';
  auxData.numProcessors = 6;
  auxData.groupSize = 1000;
  auxData.smallestSetOfCollected = 1000;
  auxData.numBACtoConsider = 455055
  auxData.thresholdForCollectingBAC = 1e-3;
  auxData.firstFlag = 1;
  auxData.createReadsFlag = 1;
  auxDataFileName = ['~/CS/tmp/structFor_',auxData.tmpFileName];
  save(auxDataFileName,'auxData')
  distributeBAC_general(auxDataFileName)

  
  clear
  name = '~/CS/BAC/dataForSim/simSec_Nbacmix_100_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads';
  origName = name
  load(name)
  
  a = find(name=='/');
  name(1:a(end)) = [];
  name(find(name=='_')) = '-';
  load ~/CS/BAC/testNew
  found
  for i=1:length(found)
    if ~isempty(found{i})
      figure(i)
      plot(found{i},correctWeight,'.');
      title([name,num2str(i)])
    end
  end
  
  
  clear
  auxData = struct;
  auxData.randomizeSeedAfterCreateReadsFlag = 1;
  auxData.reads_data_fileName = '~/CS/BAC/dataForSim/simBACdata.mat';
  auxData.readLength = 50;
  auxData.tmpFileName = 'testNewOU';
  auxData.saveName = '~/CS/BAC/testNewOU';
  auxData.numProcessors = 6;
  auxData.firstFlag = 0;
  auxData.groupSize = 1000;
  auxData.smallestSetOfCollected = 1000;
  auxData.numBACtoConsider = 455055;
  auxData.thresholdForCollectingBAC = 1e-3;
  auxData.randomizeSeedAfterCreateReadsFlag = 1;
  auxDataFileName = ['~/CS/tmp/structFor_',auxData.tmpFileName];
  save(auxDataFileName,'auxData')
  distributeBAC_general(auxDataFileName)

  
  
  % Br
  clear
  auxData = struct;
  auxData.randomizeSeedAfterCreateReadsFlag = 1;
  auxData.brFlag = 1;
  auxData.queueName = 'hour';
  auxData.reads_data_fileName = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/simBACdata.mat';
  auxData.readLength = 50;
  auxData.tmpFileName = 'testNew';
  auxData.saveName = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/testNew';
  auxData.numProcessors = 100;
  auxData.firstFlag = 0;
  auxData.groupSize = 1000;
  
  auxData.smallestSetOfCollected = 1000;
  auxData.numBACtoConsider = 455055;
  auxData.thresholdForCollectingBAC = 1e-3;
  auxDataFileName = ['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/structFor_',auxData.tmpFileName];
  save(auxDataFileName,'auxData');
  fprintf('/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/cvx/run_distributeBAC_general.sh  /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b %s\n',auxDataFileName);
  distributeBAC_general(auxDataFileName)
  
  %%%%%%%%%%%%%%%%simSec
  clear
  load(['~/CS/BAC/cellsimbac2'])
  Nbac_in_mixture = 100;
  Nreads = 2e6;
  readlen = 50;
  npower = 0.5;
  bac_dist_flag = 1;
  outfilename = '~/CS/BAC/dataForSim/simSec';
  N = 455055;
  createReadsFlag = 0;
  [strinfilename] = fun_prep_data_for_simulations_loadSpecificBAC2([],[],indsimbac_vec,Nbac_in_mixture,Nreads,readlen,npower,bac_dist_flag,[],outfilename,N,createReadsFlag);
  
  clear
  auxData = struct;
  auxData.createReadsFlag = 1;
  auxData.randomizeSeedAfterCreateReadsFlag = 1;
  auxData.brFlag = 1;
  auxData.queueName = 'hour';
  auxData.reads_data_fileName = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/simSec_Nbacmix_100_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads.mat';
  auxData.readLength = 50;
  auxData.tmpFileName = 'testSimSec';
  auxData.saveName = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/testSimSec';
  auxData.numProcessors = 100;
  auxData.firstFlag = 0;
  auxData.groupSize = 1000;
  
  auxData.smallestSetOfCollected = 1000;
  auxData.numBACtoConsider = 455055;
  auxData.thresholdForCollectingBAC = 1e-3;
  auxDataFileName = ['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/structFor_',auxData.tmpFileName];
  save(auxDataFileName,'auxData');
  fprintf('/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/cvx/run_distributeBAC_general.sh  /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b %s\n',auxDataFileName);
  distributeBAC_general(auxDataFileName)
  
  plotRes('~/CS/BAC/dataForSim/simSec_Nbacmix_100_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads.mat','~/CS/BAC/testSimSec')
  
  clear
  print('-dpdf',origName)
  
  %%%%%%% 500 close
  clear
  load(['~/CS/BAC/cellsimbac2'])
  Nbac_in_mixture = 500;
  Nreads = 2e6;
  readlen = 50;
  npower = 0.5;
  bac_dist_flag = 1;
  outfilename = '~/CS/BAC/dataForSim/simSec500';
  N = 455055;
  createReadsFlag = 0;
  [strinfilename] = fun_prep_data_for_simulations_loadSpecificBAC2([],[],indsimbac_vec,Nbac_in_mixture,Nreads,readlen,npower,bac_dist_flag,[],outfilename,N,createReadsFlag);
  
  % unix(['scp ',strinfilename,' orzuk@tin.broadinstitute.org:~/']);
  
  clear
  auxData = struct;
  auxData.createReadsFlag = 1;
  auxData.randomizeSeedAfterCreateReadsFlag = 1;
  auxData.brFlag = 1;
  auxData.queueName = 'hour';
  auxData.reads_data_fileName = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/simSec500_Nbacmix_500_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads.mat';
  auxData.readLength = 50;
  auxData.tmpFileName = 'testSimSec500';
  auxData.saveName = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/testSimSec500';
  auxData.numProcessors = 100;
  auxData.firstFlag = 0;
  auxData.groupSize = 1000;
  
  auxData.smallestSetOfCollected = 1000;
  auxData.numBACtoConsider = 455055;
  auxData.thresholdForCollectingBAC = 1e-3;
  auxDataFileName = ['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/structFor_',auxData.tmpFileName];
  save(auxDataFileName,'auxData');
  fprintf('/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/cvx/run_distributeBAC_general.sh  /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b %s\n',auxDataFileName);
  distributeBAC_general(auxDataFileName)
  
  plotRes('~/CS/BAC/dataForSim/simSec500_Nbacmix_500_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads.mat','~/CS/BAC/testSimSec500')
  
  %%%%%%%end 500 close
  
  %%%%%%%%%
  clear
  load(['~/CS/BAC/cellsimbac2'])
  Nbac_in_mixture = 100;
  Nreads = 2e7;
  readlen = 50;
  npower = 0.5;
  bac_dist_flag = 1;
  outfilename = '~/CS/BAC/dataForSim/simSec2e7';
  N = 455055;
  createReadsFlag = 0;
  [strinfilename] = fun_prep_data_for_simulations_loadSpecificBAC2([],[],indsimbac_vec,Nbac_in_mixture,Nreads,readlen,npower,bac_dist_flag,[],outfilename,N,createReadsFlag);
  
  clear
  auxData = struct;
  auxData.createReadsFlag = 1;
  auxData.randomizeSeedAfterCreateReadsFlag = 1;
  auxData.brFlag = 1;
  auxData.queueName = 'priority';
  auxData.reads_data_fileName = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/simSec2e7_Nbacmix_100_Nread_20000000_Readlen_50_Npower_05_bacdistflag_1_NoReads.mat';
  auxData.readLength = 50;
  auxData.tmpFileName = 'testSimSec2e7';
  auxData.saveName = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/testSimSec2e7';
  auxData.numProcessors = 50;
  auxData.firstFlag = 0;
  auxData.groupSize = 1000;
  auxData.smallestSetOfCollected = 1000;
  auxData.numBACtoConsider = 455055;
  auxData.thresholdForCollectingBAC = 1e-3;
  auxDataFileName = ['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/structFor_',auxData.tmpFileName];
  save(auxDataFileName,'auxData');
  fprintf('/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/cvx/run_distributeBAC_general.sh  /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b %s\n',auxDataFileName);
  distributeBAC_general(auxDataFileName)
  
  plotRes('~/CS/BAC/dataForSim/simSec2e7_Nbacmix_100_Nread_20000000_Readlen_50_Npower_05_bacdistflag_1_NoReads','~/CS/BAC/testSimSec2e7')
  
  
  %%%%%%%%%%%%%%%%end simSec

  
  
end

