% effect of 10 times
if 1==2
  clear
 auxData = struct;
 
 %userDir = getuserdir;
 userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
 auxData.repeatRandomGroups = 10; %!!!!
    %%%
    auxData.moveDependentToNextStageFlag = 1;
    auxData.tmpFileName = ['testMinFreq200MoveDependentInfiniteGroup500Repeat10Times'];
    auxData.randomizeSeedAfterCreateReadsFlag = 1;
    auxData.saveName = [userDir,'/CS/BAC/',auxData.tmpFileName];
    auxData.inifiniteNumberOfReadsFlag = 1;
    %%
    
    auxData.currDir = [userDir,'/CS/BAC/tmpRuns/',auxData.tmpFileName];
    unix(['mkdir ',auxData.currDir]);
    unix(['cd ',auxData.currDir]); 
    auxData.createReadsFlag = 1;
    auxData.brFlag = 1;
    auxData.queueName = 'hour';
    auxData.reads_data_fileName = [userDir,'/CS/BAC/dataForSim/simChangeMinFreq_Nbacmix_200_Nread_2000000_Readlen_50_Npower_01_bacdistflag_1_NoReads.mat'];
    auxData.readLength = 50;
    auxData.numProcessors = 100;
    auxData.firstFlag = 0;
    auxData.groupSize = 500;
    auxData.smallestSetOfCollected = 1000;
    auxData.numBACtoConsider = 455055;
    auxData.thresholdForCollectingBAC = 1e-3;
    auxDataFileName = [userDir,'/CS/tmp/structFor_',auxData.tmpFileName];
    save(auxDataFileName,'auxData');
    
    w1 = [sprintf([';cd ',auxData.currDir,';','/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/cvx/run_distributeBAC_general3.sh  /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b %s;'],auxDataFileName)]
    
    

end % mix of 200 in groups of 500 - do 10 time when the group is less than 20000

% 10 times infinite
if 1==2
  clear
 auxData = struct;
 
 %userDir = getuserdir;
 userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
 auxData.repeatRandomGroups = 10; %!!!!
    %%%
    auxData.moveDependentToNextStageFlag = 1;
    auxData.tmpFileName = ['testMinFreq200MoveDependentInfiniteGroup500Repeat10TimesGroupOf500'];
    auxData.randomizeSeedAfterCreateReadsFlag = 1;
    auxData.saveName = [userDir,'/CS/BAC/',auxData.tmpFileName];
    auxData.inifiniteNumberOfReadsFlag = 1;
    %%
    
    auxData.currDir = [userDir,'/CS/BAC/tmpRuns/',auxData.tmpFileName];
    unix(['mkdir ',auxData.currDir]);
    unix(['cd ',auxData.currDir]); 
    auxData.createReadsFlag = 1;
    auxData.brFlag = 1;
    auxData.queueName = 'hour';
    auxData.reads_data_fileName = [userDir,'/CS/BAC/dataForSim/simChangeMinFreq_Nbacmix_200_Nread_2000000_Readlen_50_Npower_01_bacdistflag_1_NoReads.mat'];
    auxData.readLength = 50;
    auxData.numProcessors = 100;
    auxData.firstFlag = 0;
    auxData.groupSize = 500;
    auxData.smallestSetOfCollected = 1000;
    auxData.numBACtoConsider = 455055;
    auxData.thresholdForCollectingBAC = 1e-3;
    auxDataFileName = [userDir,'/CS/tmp/structFor_',auxData.tmpFileName];
    save(auxDataFileName,'auxData');
    
    w1 = [sprintf([';cd ',auxData.currDir,';','/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/cvx/run_distributeBAC_general3.sh  /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b %s;'],auxDataFileName)]
    
    

end % mix of 200 in groups of 500 - do 10 time when the group is less than 20000

% do 10 times finite
if 1==2 
  clear
 auxData = struct;
 
 %userDir = getuserdir;
 userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
 auxData.repeatRandomGroups = 10; %!!!!
                                  %%%
 auxData.moveDependentToNextStageFlag = 1;
 auxData.tmpFileName = ['testMinFreq200MoveDependentInfiniteGroup200Repeat10TimesGroupOf500Finite'];
 auxData.randomizeSeedAfterCreateReadsFlag = 1;
 auxData.saveName = [userDir,'/CS/BAC/',auxData.tmpFileName];
 auxData.inifiniteNumberOfReadsFlag = 0;
 %%
    
 auxData.currDir = [userDir,'/CS/BAC/tmpRuns/',auxData.tmpFileName];
 unix(['mkdir ',auxData.currDir]);
 unix(['cd ',auxData.currDir]); 
 auxData.createReadsFlag = 1;
 auxData.brFlag = 1;
 auxData.queueName = 'hour';
 auxData.reads_data_fileName = [userDir,'/CS/BAC/dataForSim/simChangeMinFreq_Nbacmix_200_Nread_2000000_Readlen_50_Npower_01_bacdistflag_1_NoReads.mat'];
 auxData.readLength = 50;
 auxData.numProcessors = 100;
 auxData.firstFlag = 0;
 auxData.groupSize = 500;
 auxData.smallestSetOfCollected = 1000;
 auxData.numBACtoConsider = 455055;
 auxData.thresholdForCollectingBAC = 1e-3;
 auxDataFileName = [userDir,'/CS/tmp/structFor_',auxData.tmpFileName];
 save(auxDataFileName,'auxData');
 
 w1 = [sprintf([';cd ',auxData.currDir,';','/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/cvx/run_distributeBAC_general3.sh  /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b %s;'],auxDataFileName)]
      

end %  do 10 time finite


%%%%%%%%%%%%%%%%55

%%%%%%%5
% check preprocessing as a function of the number of reads
if 1==2
  auxData = struct;
    
    %%%
    auxData.preprocessByExistanceFlag = 1;
    auxData.moveDependentToNextStageFlag = 1;
    auxData.tmpFileName = ['testMinFreq200PreprocessingInfinite'];
    auxData.randomizeSeedAfterCreateReadsFlag = 1;
    auxData.saveName = ['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/',auxData.tmpFileName];
    auxData.inifiniteNumberOfReadsFlag = 1;
    %%
    
    auxData.currDir = ['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/tmpRuns/',auxData.tmpFileName];
    unix(['mkdir ',auxData.currDir]);
    unix(['cd ',auxData.currDir]); 
    auxData.createReadsFlag = 1;
    auxData.brFlag = 1;
    auxData.queueName = 'hour';
    auxData.reads_data_fileName = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/simChangeMinFreq_Nbacmix_200_Nread_2000000_Readlen_50_Npower_01_bacdistflag_1_NoReads.mat';
    auxData.readLength = 50;
    auxData.numProcessors = 100;
    auxData.firstFlag = 0;
    auxData.groupSize = 500;
    auxData.smallestSetOfCollected = 1000;
    auxData.numBACtoConsider = 455055;
    auxData.thresholdForCollectingBAC = 1e-3;
    auxDataFileName = ['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/structFor_',auxData.tmpFileName];
    save(auxDataFileName,'auxData');
    
    w1 = [sprintf([';cd ',auxData.currDir,';','/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/cvx/run_distributeBAC_generalOneRound.sh  /broad/software/nonfree/Linux/' ...
                    'redhat_5_x86_64/pkgs/matlab_2010b %s;'],auxDataFileName)]

    

    %%%%%%%%%%%%
    % finite
    %clear
    %userDir = getuserdir;
    userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
    load([userDir,'/CS/BAC/cellsimbac2'])
    Nbac_in_mixture = 200;
    readlen = 50;
    npower = 0.1;
    bac_dist_flag = 1;
    outfilename = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/simChangeMinFreq';
    %outfilename = '~/CS/BAC/dataForSim/simChangeMinFreq';
    N = 455055;
    createReadsFlag = 0;
    for Nreads = [10^4 10^5 10^6]
      [strinfilename] = fun_prep_data_for_simulations_loadSpecificBAC2([],[],indsimbac_vec,Nbac_in_mixture,Nreads,readlen,npower,bac_dist_flag,[],outfilename,N,createReadsFlag);
      
      
      auxData = struct;
    
      %%%
      auxData.preprocessByExistanceFlag = 1;
      auxData.moveDependentToNextStageFlag = 1;
      auxData.tmpFileName = ['testMinFreq200Preprocessing',num2str(Nreads)];
      auxData.randomizeSeedAfterCreateReadsFlag = 1;
      auxData.saveName = ['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/',auxData.tmpFileName];
      auxData.inifiniteNumberOfReadsFlag = 0;
      %%
    
      auxData.currDir = ['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/tmpRuns/',auxData.tmpFileName];
      %unix(['mkdir ',auxData.currDir]);
      %unix(['cd ',auxData.currDir]); 
      auxData.createReadsFlag = 1;
      auxData.brFlag = 1;
      auxData.queueName = 'hour';
      auxData.reads_data_fileName = strinfilename;
      auxData.readLength = 50;
      auxData.numProcessors = 100;
      auxData.firstFlag = 0;
      auxData.groupSize = 500;
      auxData.smallestSetOfCollected = 1000;
      auxData.numBACtoConsider = 455055;
      auxData.thresholdForCollectingBAC = 1e-3;
      auxDataFileName = ['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/structFor_',auxData.tmpFileName];
      
      save(auxDataFileName,'auxData');
      distributeBAC_generalOneRound(auxDataFileName);
      
      
      %w1 = [w1,sprintf([';cd ',auxData.currDir,';','/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/cvx/run_distributeBAC_generalOneRound.sh  /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b %s;'],auxDataFileName)]
    end
    
    
end % check preprocessing as a function of the number of reads

%%%%%%%%%%%%%%%%%%%%%%%%%%

% effect of 10 times
clear
plotRes('~/CS/BAC/testMinFreq200MoveDependentInfiniteGroup500Repeat10Times.mat','~/CS/BAC/dataForSim/simChangeMinFreq_Nbacmix_200_Nread_2000000_Readlen_50_Npower_01_bacdistflag_1_NoReads.mat')
print -dpdf ~/CS/BAC/res200Close_takeConcensusOfGroupsToAvoidFalseNegatives

load ~/CS/BAC/testMinFreq200MoveDependentInfiniteGroup500Repeat10Times.mat
load ~/CS/BAC/dataForSim/simChangeMinFreq_Nbacmix_200_Nread_2000000_Readlen_50_Npower_01_bacdistflag_1_NoReads.mat

a = find(abs(found{4})<10^-5 & correctWeight'>0);
% the is one which we miss  -


%%%% effect of 10 times for finite
plotRes('~/CS/BAC/testMinFreq200MoveDependentInfiniteGroup200Repeat10TimesGroupOf500Finite.mat','~/CS/BAC/dataForSim/simChangeMinFreq_Nbacmix_200_Nread_2000000_Readlen_50_Npower_01_bacdistflag_1_NoReads.mat')
print -dpdf ~/CS/BAC/res200Close_takeConcensusOfGroupsToAvoidFalseNegativesFinite

load ~/CS/BAC/testMinFreq200MoveDependentInfiniteGroup200Repeat10TimesGroupOf500Finite
load ~/CS/BAC/dataForSim/simChangeMinFreq_Nbacmix_200_Nread_2000000_Readlen_50_Npower_01_bacdistflag_1_NoReads.mat

a = find(abs(found{4})<10^-5 & correctWeight'>0);



%%%%%%%%%%%%%%
% preprocessing

clear
load  ~/CS/BAC/dataForSim/simChangeMinFreq_Nbacmix_200_Nread_2000000_Readlen_50_Npower_01_bacdistflag_1_NoReads.mat


clf
k = 1;
for Nreads = [10^4 10^5 10^6 inf]
  subplot(2,2,k)
  if isinf(Nreads)
    load ~/CS/BAC/testMinFreq200PreprocessingInfinite
  else
    load(['~/CS/BAC/testMinFreq200Preprocessing',num2str(Nreads)])
  end
  
  plot(currX,'.')
  hold on
  plot(ind_bac_in_mix,currX(ind_bac_in_mix),'ro')
  title(Nreads)
  clear currX
  k =k+1;
end
print -dpdf ~/CS/BAC/res200Close_preprocessing






