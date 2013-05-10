function crOr(basicPrefix,outFileDirName,numProcessors,bdf,nb,np,numIter,nr)
disp('works with non-ACGT data')
%userDir = getuserdir;
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
load([userDir,'/CS/BAC/cellsimbac2'])

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
      for readlen=[50]
        for ii=1:numIter
          seed = sum(100*clock);
          for Nreads=nr
            
            basicName = [basicPrefix,'_bac_dist_flag','_',num2str(bac_dist_flag),...
                         '_Nbac_in_mixture','_',num2str(Nbac_in_mixture),...
                         '_npower','_',num2str(npower*10),...
                         '_readlen','_',num2str(readlen),...
                         '_numIter','_',num2str(ii),...
                         '_Nreads','_',num2str(Nreads)];
            outfilename = [userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',basicName];
            
            [strinfilename{k}] = fun_prep_data_for_simulations_loadSpecificBAC2Or([],[],indsimbac_vec,Nbac_in_mixture,Nreads,readlen,npower,bac_dist_flag,seed,outfilename,N,createReadsFlag,1);
            saveSeed(k) = seed;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            auxData = struct;
            
            %%%
            auxData.repeatRandomGroups = numRepeatRandomGroups;
            auxData.moveDependentToNextStageFlag = 1;
          
            auxData.tmpFileName = basicName;
            auxData.randomizeSeedAfterCreateReadsFlag = 1;
            auxData.saveName = [userDir,'/CS/BAC/',outFileDirName,'/',auxData.tmpFileName];
            
            if isinf(Nreads)
              auxData.inifiniteNumberOfReadsFlag = 1;
            else
              auxData.inifiniteNumberOfReadsFlag = 0;
            end
            
            %%
            
            auxData.currDir = [userDir,'/CS/BAC/tmpRuns/',outFileDirName,'/',auxData.tmpFileName];
            unix(['mkdir ',auxData.currDir]);
            
            auxData.basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];
            auxData.basicSeqKey= [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
            
            auxData.createReadsFlag = 1;
            auxData.brFlag = 1;
            auxData.queueName = 'hour';
            auxData.reads_data_fileName = strinfilename{k};
            auxData.readLength = readlen;
            auxData.numProcessors = numProcessors;
            auxData.firstFlag = 0;
            auxData.groupSize = 1000;
            auxData.smallestSetOfCollected = 1000;
            auxData.numBACtoConsider = N;
            auxData.thresholdForCollectingBAC = 1e-3;
            auxDataFileName = [userDir,'/CS/tmp/',outFileDirName,'/','structFor_',auxData.tmpFileName];
            save(auxDataFileName,'auxData');
            
            
            
            
            
            auxDataFileName = [userDir,'/CS/tmp/',outFileDirName,'/','structFor_',basicName];
            currDir = [userDir,'/CS/BAC/tmpRuns/',outFileDirName,'/',basicName]; 
            
            currFileName = [userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','run_',basicName]
            
            curr_fid = fopen(currFileName,'w');
            fprintf(curr_fid,'#!/bin/bash\n');
            fprintf(curr_fid,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
            fprintf(curr_fid,'mkdir $TMPDIR\n');
            fprintf(curr_fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
            
            w1 = [sprintf(['cd ',currDir,';'])];
            fprintf(curr_fid,'%s\n',w1);
            w1 = ['bsub -q week -M 33554432',...
                  ' -o ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',basicName,'.o ',...
                  ' -e ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',basicName,'.e ',...
                  '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/cvx/run_distributeBAC_generalOr.sh  /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b ',auxDataFileName];
            
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
            fprintf(fid,'%s\n',currFileName);
            
            
            k = k+1;
          end
        end
      end
    end
  end
end

unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','listOfR_',basicPrefix,'*']);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','run_*'])
save([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',basicPrefix,'_list'],'totalListOfNames','totalListOfRuns')
  

