% 

% close - distant

% 0.5 uniform

% 

% create list of cases: 50: 

clear
createReadsFlag = 0;
N = 455055;
numIter = 2;
%userDir = getuserdir;
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
load([userDir,'/CS/BAC/cellsimbac2'])

fclose('all')
k = 1;
for bac_dist_flag=[0,1]
  for Nbac_in_mixture = [10,500]
    for npower=[0.5 0]
      for readlen=[50]
        for ii=1:numIter
          seed = sum(100*clock);
          for Nreads=[10^4 10^5 10^6 Inf]
            
            basicName = ['s_','bac_dist_flag','_',num2str(bac_dist_flag),...
            '_Nbac_in_mixture','_',num2str(Nbac_in_mixture),...
            '_npower','_',num2str(npower*10),...
            '_readlen','_',num2str(readlen),...
            '_numIter','_',num2str(ii),...
            '_Nreads','_',num2str(Nreads)];
            
            if k>0
              outfilename = [userDir,'/CS/BAC/dataForSim/set/',basicName];
              
              [strinfilename{k}] = fun_prep_data_for_simulations_loadSpecificBAC2([],[],indsimbac_vec,Nbac_in_mixture,Nreads,readlen,npower,bac_dist_flag,[],outfilename,N,createReadsFlag);
              saveSeed(k) = seed;
              
              
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              auxData = struct;
              
              %%%
              auxData.moveDependentToNextStageFlag = 1;
              auxData.tmpFileName = basicName;
              auxData.randomizeSeedAfterCreateReadsFlag = 1;
              auxData.saveName = [userDir,'/CS/BAC/',auxData.tmpFileName];
              
              if isinf(Nreads)
                auxData.inifiniteNumberOfReadsFlag = 1;
              else
                auxData.inifiniteNumberOfReadsFlag = 0;
              end
              
              %%
              
              auxData.currDir = [userDir,'/CS/BAC/tmpRuns/',auxData.tmpFileName];
              %unix(['mkdir ',auxData.currDir]);
              %unix(['cd ',auxData.currDir]); 
              auxData.createReadsFlag = 1;
              auxData.brFlag = 1;
              auxData.queueName = 'hour';
              auxData.reads_data_fileName = strinfilename{k};
              auxData.readLength = readlen;
              auxData.numProcessors = 100;
              auxData.firstFlag = 0;
              auxData.groupSize = 1000;
              auxData.smallestSetOfCollected = 1000;
              auxData.numBACtoConsider = 455055;
              auxData.thresholdForCollectingBAC = 1e-3;
              auxDataFileName = [userDir,'/CS/tmp/structFor_',auxData.tmpFileName];
              save(auxDataFileName,'auxData');
              
              w1 = [sprintf(['cd ',auxData.currDir,';','/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/cvx/run_distributeBAC_general3.sh  /broad/software/nonfree/Linux/' ...
                             'redhat_5_x86_64/pkgs/matlab_2010b %s;'],auxDataFileName)];
              
              
              if mod(k,10)==1
                if exist('fid')
                  fprintf(fid,'rm -fr $TMPDIR\n');
                  fclose(fid);
                  clear fid
                end
                
                fid = fopen([userDir,'/CS/BAC/dataForSim/set/listOfR_',num2str(k)],'w');
                fprintf(fid,'#!/bin/bash\n');
                fprintf(fid,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
                fprintf(fid,'mkdir $TMPDIR\n');
                fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
                fprintf(fid,'\n');
                

            
              end
              fprintf(fid,'%s\n',w1);
              
            end
            
            
            
            
            
            
            
            
            
            
            
            k = k+1;
            
            
            
          end
        end
      end
    end
  end
end
fclose(fid);

unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/set/listOfR_*']);



%%%%%%%%%%%
% create
clear
createReadsFlag = 0;
N = 455055;
numIter = 2;
%userDir = getuserdir;
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
load([userDir,'/CS/BAC/cellsimbac2'])

10 times
iter = 10
[10 200 ]

totalListOfRuns = [];
fclose('all')
k = 1;
for bac_dist_flag=[0,1]
  for Nbac_in_mixture = [10,500]
    for npower=[0.5 0]
      for readlen=[50]
        for ii=1:numIter
          seed = sum(100*clock);
          for Nreads=[10^4 10^5 10^6 Inf]
            
            basicName = ['s_','bac_dist_flag','_',num2str(bac_dist_flag),...
                         '_Nbac_in_mixture','_',num2str(Nbac_in_mixture),...
                         '_npower','_',num2str(npower*10),...
                         '_readlen','_',num2str(readlen),...
                         '_numIter','_',num2str(ii),...
                         '_Nreads','_',num2str(Nreads)];
            auxDataFileName = [userDir,'/CS/tmp/structFor_',basicName];
            currDir = [userDir,'/CS/BAC/tmpRuns/',basicName]; 
            
            currFileName = [userDir,'/CS/BAC/dataForSim/set/run_',basicName]
            
            curr_fid = fopen(currFileName,'w');
            fprintf(curr_fid,'#!/bin/bash\n');
            fprintf(curr_fid,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
            fprintf(curr_fid,'mkdir $TMPDIR\n');
            fprintf(curr_fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
            
            w1 = [sprintf(['cd ',currDir,';'])];
            fprintf(curr_fid,'%s\n',w1);
            w1 = ['bsub -q week -M 33554432',...
                  ' -o ',userDir,'/CS/BAC/dataForSim/set/',basicName,'.o ',...
                  ' -e ',userDir,'/CS/BAC/dataForSim/set/',basicName,'.e ',...
                  '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/cvx/run_distributeBAC_general3.sh  /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b ',auxDataFileName];
            
            fprintf(curr_fid,'%s\n',w1);
            fprintf(curr_fid,'rm -fr $TMPDIR\n');
            fclose(curr_fid);
            clear curr_fid
           
            if mod(k,10)==1
              if exist('fid')
                fprintf(fid,'rm -fr $TMPDIR\n');
                fclose(fid);
                clear fid
              end
              
              fid = fopen([userDir,'/CS/BAC/dataForSim/set/listOfR_',num2str(k)],'w');
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


unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/set/listOfR_*']);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/set/run_*'])


% 10 per CPU
% check which has not ran
% load the results and compare 


%%%%%%%%%%%%%%%%%%%%%

plotRes('~/CS/BAC/s_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_10000','~/CS/BAC/dataForSim/set/s_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_10000_Nbacmix_10_Nread_10000_Readlen_50_Npower_05_bacdistflag_0_NoReads')

plotRes('~/CS/BAC/s_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_100000','~/CS/BAC/dataForSim/set/s_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_100000_Nbacmix_10_Nread_100000_Readlen_50_Npower_05_bacdistflag_0_NoReads')

plotRes('~/CS/BAC/s_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_1000000','~/CS/BAC/dataForSim/set/s_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_1000000_Nbacmix_10_Nread_1000000_Readlen_50_Npower_05_bacdistflag_0_NoReads')




ran 1
11
the first in 21
31: in upper left window
21,51,41 runs in bsub set/log_listOfR21 51,41
61: runs in bsub set/log_listOfR61

new
1
11 - done
21 - running at 15:30 
31 - running
41 - running
51 - runng
61 - runnong

41
51
61

running - started at 6am 
1 - done
11
41
61


cd /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/set; bsub -q week -M 33554432 -o /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/set/log_listR11.o -e /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/set/log_listR11.e /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/set/listOfR_11

cd /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/set; bsub -q week -M 33554432 -o /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/set/log_listR21.o -e /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/set/log_listR21.e /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/set/listOfR_21

cd /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/set; bsub -q week -M 33554432 -o /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/set/log_listR31.o -e /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/set/log_listR31.e /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/set/listOfR_31

cd /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/set; bsub -q week -M 33554432 -o /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/set/log_listR41.o -e /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/set/log_listR41.e /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/set/listOfR_41

cd /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/set; bsub -q week -M 33554432 -o /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/set/log_listR51.o -e /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/set/log_listR51.e /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/set/listOfR_51

cd /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/set; bsub -q week -M 33554432 -o /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/set/log_listR61.o -e /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/set/log_listR61.e /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/set/listOfR_61


s_bac_dist_flag_0_Nbac_in_mixture_500_npower_0_readlen_50_numIter_2_Nreads_100000_Nbacmix_500_Nread_100000_Readlen_50_Npower_0_bacdistflag_0_NoReads

cd /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC;scp s_bac_dist_flag* shental@sol.cslab.openu.ac.il:~/CS/BAC/



s_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_2_Nreads_1000000

distributeBAC_general3('/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/structFor_s_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_10000')



/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/tmpRuns/s_bac_dist_flag_0_Nbac_in_mixture_10_npower_0_readlen_50_numIter_1_Nreads_Inf/run_s_bac_dist_flag_0_Nbac_in_mixture_10_npower_0_readlen_50_numIter_1_Nreads_Inf_k_1_logOfUnsuccessful3

/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/tmpRuns/s_bac_dist_flag_0_Nbac_in_mixture_10_npower_0_readlen_50_numIter_2_Nreads_10000/run_s_bac_dist_flag_0_Nbac_in_mixture_10_npower_0_readlen_50_numIter_2_Nreads_10000_k_1_logOfUnsuccessful3.t-xt



/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/tmpRuns/s_bac_dist_flag_0_Nbac_in_mixture_10_npower_0_readlen_50_numIter_1_Nreads_Inf/run_s_bac_dist_flag_0_Nbac_in_mixture_10_npower_0_readlen_50_numIter_1_Nreads_Inf_k_1_logOfUnsuccessful3.txt


d = str2mat('s_bac_dist_flag_0_Nbac_in_mixture_500_npower_5_readlen_50_numIter_1_Nreads_Inf',...
's_bac_dist_flag_0_Nbac_in_mixture_10_npower_0_readlen_50_numIter_2_Nreads_10000',...
's_bac_dist_flag_0_Nbac_in_mixture_500_npower_5_readlen_50_numIter_1_Nreads_10000',...
's_bac_dist_flag_0_Nbac_in_mixture_500_npower_5_readlen_50_numIter_1_Nreads_100000',...
's_bac_dist_flag_0_Nbac_in_mixture_10_npower_0_readlen_50_numIter_2_Nreads_1000000',...
's_bac_dist_flag_0_Nbac_in_mixture_10_npower_0_readlen_50_numIter_1_Nreads_Inf',...
's_bac_dist_flag_0_Nbac_in_mixture_10_npower_0_readlen_50_numIter_2_Nreads_Inf',...
's_bac_dist_flag_0_Nbac_in_mixture_10_npower_0_readlen_50_numIter_2_Nreads_100000',...
's_bac_dist_flag_1_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_10000',...
's_bac_dist_flag_1_Nbac_in_mixture_10_npower_0_readlen_50_numIter_1_Nreads_Inf',...
's_bac_dist_flag_1_Nbac_in_mixture_10_npower_5_readlen_50_numIter_2_Nreads_Inf',...
's_bac_dist_flag_1_Nbac_in_mixture_10_npower_5_readlen_50_numIter_2_Nreads_10000',...
's_bac_dist_flag_1_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_1000000',...
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_0_readlen_50_numIter_2_Nreads_10000',...
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_0_readlen_50_numIter_2_Nreads_100000',...
's_bac_dist_flag_0_Nbac_in_mixture_500_npower_5_readlen_50_numIter_2_Nreads_1000000',...
's_bac_dist_flag_1_Nbac_in_mixture_10_npower_0_readlen_50_numIter_2_Nreads_10000',...
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_0_readlen_50_numIter_2_Nreads_1000000',...
's_bac_dist_flag_0_Nbac_in_mixture_500_npower_5_readlen_50_numIter_2_Nreads_Inf',...
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_5_readlen_50_numIter_2_Nreads_10000',...
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_0_readlen_50_numIter_1_Nreads_100000',...
's_bac_dist_flag_1_Nbac_in_mixture_10_npower_0_readlen_50_numIter_2_Nreads_100000',...
's_bac_dist_flag_0_Nbac_in_mixture_10_npower_0_readlen_50_numIter_1_Nreads_1000000',...
's_bac_dist_flag_0_Nbac_in_mixture_500_npower_0_readlen_50_numIter_2_Nreads_100000',...
's_bac_dist_flag_0_Nbac_in_mixture_500_npower_0_readlen_50_numIter_1_Nreads_10000',...
's_bac_dist_flag_1_Nbac_in_mixture_10_npower_0_readlen_50_numIter_2_Nreads_Inf',...
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_5_readlen_50_numIter_1_Nreads_Inf',...
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_5_readlen_50_numIter_2_Nreads_Inf',...
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_0_readlen_50_numIter_1_Nreads_Inf',...
's_bac_dist_flag_0_Nbac_in_mixture_500_npower_0_readlen_50_numIter_2_Nreads_10000',...
's_bac_dist_flag_0_Nbac_in_mixture_500_npower_0_readlen_50_numIter_2_Nreads_Inf',...
's_bac_dist_flag_1_Nbac_in_mixture_10_npower_0_readlen_50_numIter_1_Nreads_1000000',...
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_0_readlen_50_numIter_1_Nreads_1000000',...
's_bac_dist_flag_0_Nbac_in_mixture_500_npower_0_readlen_50_numIter_1_Nreads_100000',...
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_5_readlen_50_numIter_1_Nreads_1000000',...
's_bac_dist_flag_0_Nbac_in_mixture_500_npower_5_readlen_50_numIter_2_Nreads_10000',...
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_5_readlen_50_numIter_2_Nreads_1000000',...
's_bac_dist_flag_0_Nbac_in_mixture_500_npower_0_readlen_50_numIter_1_Nreads_1000000',...
's_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_100000',...
's_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_1000000',...
's_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_Inf',...
's_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_2_Nreads_10000',...
's_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_2_Nreads_100000',...
's_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_2_Nreads_1000000',...
's_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_2_Nreads_Inf',...
's_bac_dist_flag_0_Nbac_in_mixture_10_npower_0_readlen_50_numIter_1_Nreads_10000',...
's_bac_dist_flag_0_Nbac_in_mixture_10_npower_0_readlen_50_numIter_1_Nreads_100000',...
's_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_10000')

for i=1:size(d,1)
  close all
  [Nbac_in_mixture,readlen,npower,bac_dist_flag,Nreads]=findParametersBasedResFile(d(i,:));
  
  
  plotRes(['~/CS/BAC/',deblank(d(i,:)),'.mat'],['~/CS/BAC/dataForSim/set/',deblank(d(i,:)),'_Nbacmix_',num2str(Nbac_in_mixture),'_Nread_',num2str(Nreads),'_Readlen_',num2str(readlen),'_Npower_',npower,'_bacdistflag_',num2str(bac_dist_flag),'_NoReads'])

  pause
end


l = str2mat('s_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_10000',...
's_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_100000',...
's_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_1000000',...
's_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_Inf',...
's_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_2_Nreads_10000',...
's_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_2_Nreads_100000',...
's_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_2_Nreads_1000000',...
's_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_2_Nreads_Inf',...
's_bac_dist_flag_0_Nbac_in_mixture_10_npower_0_readlen_50_numIter_1_Nreads_10000',...
's_bac_dist_flag_0_Nbac_in_mixture_10_npower_0_readlen_50_numIter_1_Nreads_100000',...
's_bac_dist_flag_0_Nbac_in_mixture_500_npower_5_readlen_50_numIter_2_Nreads_10000',...
's_bac_dist_flag_0_Nbac_in_mixture_500_npower_5_readlen_50_numIter_2_Nreads_100000',...
's_bac_dist_flag_0_Nbac_in_mixture_500_npower_5_readlen_50_numIter_2_Nreads_1000000',...
's_bac_dist_flag_0_Nbac_in_mixture_500_npower_5_readlen_50_numIter_2_Nreads_Inf',...
's_bac_dist_flag_0_Nbac_in_mixture_500_npower_0_readlen_50_numIter_1_Nreads_10000',...
's_bac_dist_flag_0_Nbac_in_mixture_500_npower_0_readlen_50_numIter_1_Nreads_100000',...
's_bac_dist_flag_0_Nbac_in_mixture_500_npower_0_readlen_50_numIter_1_Nreads_1000000',...
's_bac_dist_flag_0_Nbac_in_mixture_500_npower_0_readlen_50_numIter_1_Nreads_Inf',...
's_bac_dist_flag_0_Nbac_in_mixture_500_npower_0_readlen_50_numIter_2_Nreads_10000',...
's_bac_dist_flag_0_Nbac_in_mixture_500_npower_0_readlen_50_numIter_2_Nreads_100000',...
's_bac_dist_flag_0_Nbac_in_mixture_500_npower_0_readlen_50_numIter_2_Nreads_1000000',...
's_bac_dist_flag_0_Nbac_in_mixture_500_npower_0_readlen_50_numIter_2_Nreads_Inf',...
's_bac_dist_flag_1_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_10000',...
's_bac_dist_flag_1_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_100000',... 
's_bac_dist_flag_1_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_1000000',...
's_bac_dist_flag_1_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_Inf',...
's_bac_dist_flag_1_Nbac_in_mixture_10_npower_5_readlen_50_numIter_2_Nreads_10000',...
's_bac_dist_flag_1_Nbac_in_mixture_10_npower_5_readlen_50_numIter_2_Nreads_100000',...
's_bac_dist_flag_1_Nbac_in_mixture_10_npower_5_readlen_50_numIter_2_Nreads_1000000',...
's_bac_dist_flag_1_Nbac_in_mixture_10_npower_5_readlen_50_numIter_2_Nreads_Inf',...
's_bac_dist_flag_1_Nbac_in_mixture_10_npower_0_readlen_50_numIter_1_Nreads_10000',...  
's_bac_dist_flag_1_Nbac_in_mixture_10_npower_0_readlen_50_numIter_1_Nreads_100000',...
's_bac_dist_flag_1_Nbac_in_mixture_10_npower_0_readlen_50_numIter_1_Nreads_1000000',...
's_bac_dist_flag_1_Nbac_in_mixture_10_npower_0_readlen_50_numIter_1_Nreads_Inf',...
's_bac_dist_flag_1_Nbac_in_mixture_10_npower_0_readlen_50_numIter_2_Nreads_10000',...
's_bac_dist_flag_1_Nbac_in_mixture_10_npower_0_readlen_50_numIter_2_Nreads_100000',...
's_bac_dist_flag_1_Nbac_in_mixture_10_npower_0_readlen_50_numIter_2_Nreads_1000000',...
's_bac_dist_flag_1_Nbac_in_mixture_10_npower_0_readlen_50_numIter_2_Nreads_Inf',...
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_5_readlen_50_numIter_1_Nreads_10000',...
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_5_readlen_50_numIter_1_Nreads_100000',...
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_5_readlen_50_numIter_1_Nreads_1000000',...
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_5_readlen_50_numIter_1_Nreads_Inf',... 
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_5_readlen_50_numIter_2_Nreads_10000',...
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_5_readlen_50_numIter_2_Nreads_100000',...
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_5_readlen_50_numIter_2_Nreads_1000000',...
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_5_readlen_50_numIter_2_Nreads_Inf',...  
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_0_readlen_50_numIter_1_Nreads_10000',... 
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_0_readlen_50_numIter_1_Nreads_100000',...
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_0_readlen_50_numIter_1_Nreads_1000000',...
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_0_readlen_50_numIter_1_Nreads_Inf',...
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_0_readlen_50_numIter_2_Nreads_10000',... 
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_0_readlen_50_numIter_2_Nreads_100000',... 
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_0_readlen_50_numIter_2_Nreads_1000000',...
's_bac_dist_flag_1_Nbac_in_mixture_500_npower_0_readlen_50_numIter_2_Nreads_Inf')  



setdiff(l,d,'rows')