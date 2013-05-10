function crGeneralFourth_splitCorrection(basicPrefix,outFileDirName,numProcessors,bdf,nb,np,numIter,nr,addNoiseFlagInput,correctMeasurementForReadErrorsFlagInput,readLength,slave_queueName,master_queueName,upperLimitForRepeat,repeatWhenLowerThanThisValue,batchSize,brFlag,mkdirFlag,thresholdForCollectingBAC,groupSize,dataSetPrimers750Flag,CorrectReadsNew,doL1Flag)

if length(numProcessors)>1
  disp(['numProcessors is: ',num2str(numProcessors{1}) ,' without correction; and ',num2str(numProcessors{2}),' for correction'])
  disp('refers only to ifinte reads!. crGeneralFourth.m')
end

if ~exist('doL1Flag')
  doL1Flag = 0;
end

if ~exist('dataSetPrimers750Flag')
  dataSetPrimers750Flag = 0;
end

if ~exist('thresholdForCollectingBAC')
  thresholdForCollectingBAC = 10^-3;
end

if ~exist('groupSize')
  groupSize = 1000;
end

if ~exist('justListFlag')
  justListFlag = 0;
end

disp('works with non-ACGT data')

userDirForFile = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';

if brFlag
  userDir = userDirForFile;
else
  userDir = getuserdir;
end
load([userDir,'/CS/BAC/cellsimbac2'])

mkdirW = [];

createReadsFlag = 0;
if dataSetPrimers750Flag==0
  N = 410849;
else
  N = 212040;
end

numRepeatRandomGroups = 10;

totalListOfRuns_noCorrection = [];
totalListOfNames_noCorrection = [];
totalListOfRuns_withCorrection = [];
totalListOfNames_withCorrection = [];

fclose('all')
k_noCorrection = 1;
k_withCorrection = 1;

numCPU_counter = inf;

for bac_dist_flag=bdf
  for Nbac_in_mixture = nb
    for npower=np
      for readlen=readLength
        
        for ii=1:numIter
          seed = sum(100*clock);
          for Nreads=nr
            for addNoiseFlag=addNoiseFlagInput
              tmpCorrect = correctMeasurementForReadErrorsFlagInput;
              if addNoiseFlag==0
                delCorrect = find(tmpCorrect==1);
                tmpCorrect(delCorrect) = [];
              end
              
              if addNoiseFlag==1 & isinf(Nreads)
                tmpCorrect = [];
              end
              
              for correctMeasurementForReadErrorsFlag=tmpCorrect
                basicName = [basicPrefix,'_bac_dist_flag','_',num2str(bac_dist_flag),...
                             '_Nbac_in_mixture','_',num2str(Nbac_in_mixture),...
                             '_npower','_',num2str(npower*10),...
                             '_readlen','_',num2str(readlen),...
                             '_numIter','_',num2str(ii),...
                             '_Nreads','_',num2str(Nreads),...
                             '_Noise_',num2str(addNoiseFlag),...
                             '_Correction_',num2str(correctMeasurementForReadErrorsFlag)...
                            ];
                
                
                
                
                
                outfilename = [userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',basicName];
                
                currDirForRun = [userDir,'/CS/BAC/tmpRuns/',outFileDirName,'/',basicName];
                
                
                
                if mkdirFlag==0
                  [strinfilename] = fun_prep_data_for_simulations_loadSpecificBAC2Or([],[],indsimbac_vec,Nbac_in_mixture,Nreads,readlen,npower,bac_dist_flag,seed,outfilename,N,createReadsFlag,1);
                  %saveSeed(k) = seed;
                  
                  
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  auxData = struct;
                  
                  auxData.userDir = userDir;
                  auxData.dataSetPrimers750Flag = dataSetPrimers750Flag;
                  auxData.doL1Flag = doL1Flag;
                  %%%
                  auxData.repeatRandomGroups = numRepeatRandomGroups;
                  auxData.upperLimitForRepeat = upperLimitForRepeat;
                  auxData.repeatWhenLowerThanThisValue = repeatWhenLowerThanThisValue;
                  auxData.batchSize = batchSize;
                  
                  auxData.moveDependentToNextStageFlag = 1;
                  
                  auxData.tmpFileName = basicName;
                  auxData.randomizeSeedAfterCreateReadsFlag = 1;
                  auxData.saveName = [userDirForFile,'/CS/BAC/',outFileDirName,'/',auxData.tmpFileName];
                  
                  
                  if isinf(Nreads)
                    auxData.inifiniteNumberOfReadsFlag = 1;
                  else
                    auxData.inifiniteNumberOfReadsFlag = 0;
                  end
                  
                  %%
                  
                  auxData.currDir = [userDirForFile,'/CS/BAC/tmpRuns/',outFileDirName,'/',auxData.tmpFileName];
                  
                  if dataSetPrimers750Flag==0
                    disp('full dataSet')
                    auxData.basicSeqNameDir = [userDirForFile,'/CS/BAC/datNoNonACGT/packed64/'];
                    auxData.basicSeqKey= [userDirForFile,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
                  else
                    disp('primers 750 data')
                    auxData.basicSeqNameDir = [userDirForFile,'/CS/BAC/primers750/datNoNonACGT/packed64/'];
                    auxData.basicSeqKey= [userDirForFile,'/CS/BAC/primers750/datNoNonACGT/keyNoNonACGT_primers750'];
                  end
                  
                  
                  auxData.createReadsFlag = 1;
                  auxData.brFlag = 1;
                  auxData.queueName = slave_queueName;
                  %keyboard
                  if brFlag
                    auxData.reads_data_fileName = strinfilename;
                  else
                    tmpStr = strrep(strinfilename,userDir,userDirForFile);
                    auxData.reads_data_fileName = tmpStr;
                  end
                  
                  auxData.readLength = readlen;
                  auxData.numProcessors = numProcessors{1};
                  auxData.firstFlag = 0;
                  auxData.groupSize = groupSize;
                  auxData.smallestSetOfCollected = 1000;
                  auxData.numBACtoConsider = N;
                  auxData.thresholdForCollectingBAC = thresholdForCollectingBAC;
                  auxDataFileName = [userDir,'/CS/tmp/',outFileDirName,'/','structFor_',auxData.tmpFileName];
                  
                  
                  
                  
                  if addNoiseFlag 
                    
                    ErrorStruct = []; ErrorStruct.error_model = 'exponential';
                    ErrorStruct.baseline_error = 0.005; 
                    ErrorStruct.final_error = 0.03; 
                    
                    p_one_nuc_error_mat = [ -1   0.3  0.22  0.18
                                        0.5    -1  0.22   0.6
                                        0.35  0.15    -1  0.22
                                        0.15  0.55  0.56   -1]; % , 16, length(p));
                    
                    ErrorStruct.p_one_nuc_error_mat = p_one_nuc_error_mat;
                    auxData.ErrorStruct = ErrorStruct;
                    
                    auxData.addNoiseFlag = 1;
                    
                    if correctMeasurementForReadErrorsFlag==1
                      auxData.correctMeasurementForReadErrorsFlag = 1;
                      auxData.CorrectReadsNew = CorrectReadsNew;
                      
                      % change 19.1.12
                      auxData.numProcessors = numProcessors{2}(auxData.readLength);
                      % end change 19.1.12
                    else
                      auxData.correctMeasurementForReadErrorsFlag = 0;
                    
                    end
                  else % addNoiseFlag
                    auxData.addNoiseFlag = 0;
                    auxData.correctMeasurementForReadErrorsFlag = 0;
                  end
                  
                  auxData.saveFullDataFlag = 0;
                  
                  
                  
                  
                  if (auxData.inifiniteNumberOfReadsFlag==0 & ...
                     auxData.readLength>=400 & ...
                     Nreads>10^5) | ...
                     (auxData.inifiniteNumberOfReadsFlag==0 & ...
                     Nreads>10^6)   
                    
                    auxData.createReadsAndQuit=1;
                    auxData.loadSavedReadsForLargeReadLengthAndManyReads = 1;
                    
                    %keyboard
                    save([auxDataFileName,'_createReads'],'auxData')
                    
                    auxData.createReadsAndQuit=0;
                    auxData.loadSavedReadsForLargeReadLengthAndManyReads=1;
                    
                    save([auxDataFileName],'auxData')
                    
                  else
                    save(auxDataFileName,'auxData'); 
                  end
                    
                  
                  auxDataFileName = [userDirForFile,'/CS/tmp/',outFileDirName,'/','structFor_',basicName];
                  currDir = [userDirForFile,'/CS/BAC/tmpRuns/',outFileDirName,'/',basicName]; 
                  
                  currFileName = [userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','run_',basicName];
                  currFileNameForFile = [userDirForFile,'/CS/BAC/dataForSim/',outFileDirName,'/','run_',basicName];
                  
                  curr_fid = fopen(currFileName,'w');
                  fprintf(curr_fid,'#!/bin/bash\n');
                  fprintf(curr_fid,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
                  fprintf(curr_fid,'mkdir $TMPDIR\n');
                  fprintf(curr_fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
                  
                  w1 = [sprintf(['cd ',currDir,';'])];
                  fprintf(curr_fid,'%s\n',w1);
                  w1 = ['bsub -c 1440 -q ',master_queueName,...
                        ' -o ',userDirForFile,'/CS/BAC/dataForSim/',outFileDirName,'/',basicName,'.o ',...
                        ' -e ',userDirForFile,'/CS/BAC/dataForSim/',outFileDirName,'/',basicName,'.e ',...
                        '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/cvx/run_distributeBAC_generalOrFourth.sh  /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b ',auxDataFileName];
                  % deleted ' -M 33554432',
                  
                  fprintf(curr_fid,'%s\n',w1);
                  fprintf(curr_fid,'rm -fr $TMPDIR\n');
                  fclose(curr_fid);
                  clear curr_fid
                  
                  if correctMeasurementForReadErrorsFlag==0
                    totalListOfRuns_noCorrection{k_noCorrection} = w1;
                    totalListOfNames_noCorrection{k_noCorrection} = basicName;
                    
                    if mod(k_noCorrection,10)==1
                      if exist('fid_noCorrection')
                        fprintf(fid_noCorrection,'rm -fr $TMPDIR\n');
                        fclose(fid_noCorrection);
                        clear fid_noCorrection
                      end
                      
                      fid_noCorrection = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','listOfR_noCorrection_',basicPrefix,'_',num2str(k_noCorrection)],'w');
                      fprintf(fid_noCorrection,'#!/bin/bash\n');
                      fprintf(fid_noCorrection,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
                      fprintf(fid_noCorrection,'mkdir $TMPDIR\n');
                      fprintf(fid_noCorrection,'export MCR_CACHE_ROOT=$TMPDIR\n');
                      fprintf(fid_noCorrection,'\n');
                      
                    end
                    
                    k_noCorrection = k_noCorrection+1;
                    fprintf(fid_noCorrection,'%s\n',currFileNameForFile);
                  else % with correction
                   
                    totalListOfRuns_withCorrection{k_withCorrection} = w1;
                    totalListOfNames_withCorrection{k_withCorrection} = basicName;

                    % split at 500 CPU
                    if numCPU_counter>500 %if mod(k_withCorrection,10)==1
                      numCPU_counter = 0;
                      
                      if exist('fid_withCorrection')
                        fprintf(fid_withCorrection,'sleep 4h;\n');
                        fprintf(fid_withCorrection,'rm -fr $TMPDIR\n');
                        fclose(fid_withCorrection);
                        clear fid_withCorrection
                      end
                      
                      fid_withCorrection = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','listOfR_withCorrection_',basicPrefix,'_',num2str(k_withCorrection)],'w');
                      fprintf(fid_withCorrection,'#!/bin/bash\n');
                      fprintf(fid_withCorrection,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
                      fprintf(fid_withCorrection,'mkdir $TMPDIR\n');
                      fprintf(fid_withCorrection,'export MCR_CACHE_ROOT=$TMPDIR\n');
                      fprintf(fid_withCorrection,'\n');
                      
                    end
                    
                    
                    k_withCorrection = k_withCorrection+1;
                    numCPU_counter = numCPU_counter+auxData.numProcessors;
                    fprintf(fid_withCorrection,'%s\n',currFileNameForFile);
                    
                    
                  end
                 
                  
                else
                  mkdirW = [mkdirW,';mkdir ',currDirForRun];
                  
                  
                end %mkdirFlag
                
                
                
              end % correctMeasurementForReadErrorsFlag
              
            end % addNoiseFlag
            
            
          end
        end
      end
    end
  end
end

if mkdirFlag
  %keyboard
  eval(mkdirW)
  return
end

if exist('fid_noCorrection')
  fprintf(fid_noCorrection,'rm -fr $TMPDIR\n');
  fclose(fid_noCorrection);
  clear fid_noCorrection
end

if exist('fid_withCorrection')
  fprintf(fid_withCorrection,'rm -fr $TMPDIR\n');
  fclose(fid_withCorrection);
  clear fid_withCorrection
end

unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','listOfR_*',basicPrefix,'*']);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','run_*'])
save([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',basicPrefix,'_list'],'totalListOfNames_withCorrection','totalListOfRuns_withCorrection','totalListOfNames_noCorrection','totalListOfRuns_noCorrection')
 

