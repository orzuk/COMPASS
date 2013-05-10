function crGeneralThird(basicPrefix,outFileDirName,numProcessors,bdf,nb,np,numIter,nr,addNoiseFlagInput,correctMeasurementForReadErrorsFlagInput,readLength,slave_queueName,master_queueName,upperLimitForRepeat,repeatWhenLowerThanThisValue,batchSize,brFlag,mkdirFlag)

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
N = 410849;
numRepeatRandomGroups = 10;

totalListOfRuns = [];
totalListOfNames = [];
fclose('all')
k = 1;
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
                  [strinfilename{k}] = fun_prep_data_for_simulations_loadSpecificBAC2Or([],[],indsimbac_vec,Nbac_in_mixture,Nreads,readlen,npower,bac_dist_flag,seed,outfilename,N,createReadsFlag,1);
                  saveSeed(k) = seed;
                  
                  
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  auxData = struct;
                  
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
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  auxData.basicSeqNameDir = [userDirForFile,'/CS/BAC/datNoNonACGT/packed64/'];
                  auxData.basicSeqKey= [userDirForFile,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
                  
                  auxData.createReadsFlag = 1;
                  auxData.brFlag = 1;
                  auxData.queueName = slave_queueName;
                  %keyboard
                  if brFlag
                    auxData.reads_data_fileName = strinfilename{k};
                  else
                    tmpStr = strrep(strinfilename{k},userDir,userDirForFile);
                    auxData.reads_data_fileName = tmpStr;
                  end
                  
                  auxData.readLength = readlen;
                  auxData.numProcessors = numProcessors;
                  auxData.firstFlag = 0;
                  auxData.groupSize = 1000;
                  auxData.smallestSetOfCollected = 1000;
                  auxData.numBACtoConsider = N;
                  auxData.thresholdForCollectingBAC = 1e-3;
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
                    else
                      auxData.correctMeasurementForReadErrorsFlag = 0;
                    end
                  else % addNoiseFlag
                    auxData.addNoiseFlag = 0;
                    auxData.correctMeasurementForReadErrorsFlag = 0;
                  end
                  
                  auxData.saveFullDataFlag = 0;
                  
                  save(auxDataFileName,'auxData');
                  
                  
                  
                  
                  
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
                        '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/cvx/run_distributeBAC_generalOrThird.sh  /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b ',auxDataFileName];
                  % deleted ' -M 33554432',
                  
                  fprintf(curr_fid,'%s\n',w1);
                  fprintf(curr_fid,'rm -fr $TMPDIR\n');
                  fclose(curr_fid);
                  clear curr_fid
                  
                  totalListOfRuns{k} = w1;
                  totalListOfNames{k} = basicName;
                  
                  
                  if mod(k,10)==1
                    if exist('fid')
                      fprintf(fid,'rm -fr $TMPDIR\n');
                      fclose(fid);
                      clear fid
                    end
                    
                    fid = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','listOfR_',basicPrefix,'_',num2str(k)],'w');
                    fprintf(fid,'#!/bin/bash\n');
                    fprintf(fid,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
                    fprintf(fid,'mkdir $TMPDIR\n');
                    fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
                    fprintf(fid,'\n');
                    
                  end
                  
                  
                  
                  
                  k = k+1;
                  
                  
                  
                  fprintf(fid,'%s\n',currFileNameForFile);
                  
                  
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

if exist('fid')
  fprintf(fid,'rm -fr $TMPDIR\n');
  fclose(fid);
  clear fid
end


unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','listOfR_',basicPrefix,'*']);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','run_*'])
save([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',basicPrefix,'_list'],'totalListOfNames','totalListOfRuns')
  

