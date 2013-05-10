%%%%%%%%%%%
% create
clear

%userDir = getuserdir;
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
load([userDir,'/CS/BAC/cellsimbac2'])

createReadsFlag = 0;
N = 455055;
numIter = 10;
numRepeatRandomGroups = 10;

totalListOfRuns = [];
totalListOfNames = [];
fclose('all')
k = 1;
for bac_dist_flag=[0,1]
  for Nbac_in_mixture = [10,200,500]
    for npower=[0.5 0]
      for readlen=[50]
        for ii=1:numIter
          seed = sum(100*clock);
          for Nreads=[10^4 10^5 10^6 Inf]
            
            basicName = ['ss_','bac_dist_flag','_',num2str(bac_dist_flag),...
                         '_Nbac_in_mixture','_',num2str(Nbac_in_mixture),...
                         '_npower','_',num2str(npower*10),...
                         '_readlen','_',num2str(readlen),...
                         '_numIter','_',num2str(ii),...
                         '_Nreads','_',num2str(Nreads)];
            outfilename = [userDir,'/CS/BAC/dataForSim/set/sec/',basicName];
            
            [strinfilename{k}] = fun_prep_data_for_simulations_loadSpecificBAC2([],[],indsimbac_vec,Nbac_in_mixture,Nreads,readlen,npower,bac_dist_flag,[],outfilename,N,createReadsFlag);
            saveSeed(k) = seed;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            auxData = struct;
            
            %%%
            auxData.repeatRandomGroups = numRepeatRandomGroups;
            auxData.moveDependentToNextStageFlag = 1;
          
            auxData.tmpFileName = basicName;
            auxData.randomizeSeedAfterCreateReadsFlag = 1;
            auxData.saveName = [userDir,'/CS/BAC/sec/',auxData.tmpFileName];
            
            if isinf(Nreads)
              auxData.inifiniteNumberOfReadsFlag = 1;
            else
              auxData.inifiniteNumberOfReadsFlag = 0;
            end
            
            %%
            
            auxData.currDir = [userDir,'/CS/BAC/tmpRuns/sec/',auxData.tmpFileName];
            unix(['mkdir ',auxData.currDir]);
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
            auxDataFileName = [userDir,'/CS/tmp/sec/structFor_',auxData.tmpFileName];
            save(auxDataFileName,'auxData');
            
            
            
            
            
            auxDataFileName = [userDir,'/CS/tmp/sec/structFor_',basicName];
            currDir = [userDir,'/CS/BAC/tmpRuns/sec/',basicName]; 
            
            currFileName = [userDir,'/CS/BAC/dataForSim/set/sec/run_',basicName]
            
            curr_fid = fopen(currFileName,'w');
            fprintf(curr_fid,'#!/bin/bash\n');
            fprintf(curr_fid,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
            fprintf(curr_fid,'mkdir $TMPDIR\n');
            fprintf(curr_fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
            
            w1 = [sprintf(['cd ',currDir,';'])];
            fprintf(curr_fid,'%s\n',w1);
            w1 = ['bsub -q week -M 33554432',...
                  ' -o ',userDir,'/CS/BAC/dataForSim/set/sec/',basicName,'.o ',...
                  ' -e ',userDir,'/CS/BAC/dataForSim/set/sec/',basicName,'.e ',...
                  '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/cvx/run_distributeBAC_general3.sh  /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b ',auxDataFileName];
            
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
              
              fid = fopen([userDir,'/CS/BAC/dataForSim/set/sec/listOfR_',num2str(k)],'w');
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


unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/set/sec/listOfR_*']);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/set/sec/run_*'])
save([userDir,'/CS/BAC/dataForSim/set/sec/ss_list'],'totalListOfNames','totalListOfRuns')

% run
r = 11:10:471;


s = r(1:5:end)

w = [];
for i=s
  w = [w,'./listOfR_',num2str(i),';sleep 10;'];
end

listOfR_1 - running

start 04:45
maybe runing:s

ss = setdiff(r,s);


w = ['sleep 3h;'];
for i=ss(21:30)
  w = [w,'./listOfR_',num2str(i),';sleep 10;'];
end
w = [w,'sleep 3h;']
for i=ss(31:end)
  w = [w,'./listOfR_',num2str(i),';sleep 10;'];
end



% 10 per CPU
% check which has not ran
% load the results and compare 

userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
cd([userDir,'/CS/BAC/dataForSim/set/sec'])
unix(['cd ',userDir,'/CS/BAC/dataForSim/set/sec',';ls -l *.e | grep 368 | cut -d" " -f10| cut -f1 -d"." > listOfErrors'])

fid = fopen([userDir,'/CS/BAC/dataForSim/set/sec/listOfErrors'],'r');
a = fscanf(fid,'%s\n');
fclose(fid);

fid = fopen([userDir,'/CS/BAC/dataForSim/set/sec/rerunErrors'],'w');
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
fprintf(fid,'mkdir $TMPDIR\n');
fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
fprintf(fid,'\n');


b = findstr(a,'ss_');
b = [b,length(a)+1];
for i=1:length(b)-1
  for j=1:length(totalListOfNames)
    if strcmp(a(b(i):b(i+1)-1),totalListOfNames{j})
      totalListOfRuns{j}
      fprintf(fid,'sleep 30;\n')
      fprintf(fid,'%s\n',totalListOfRuns{j});
    end
  end
end
fclose(fid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
userDir = getuserdir;
numIter = 10;
Nbac_in_mixture = [10,200,500]
Nreads=[10^4 10^5 10^6 Inf]
for bac_dist_flag=[0,1]
  for nb = 1:length(Nbac_in_mixture)
    for npower=[0.5 0]
      for readlen=[50]
        for nr = 1:length(Nreads)
          for ii=1:numIter
            
          
            if isinf(Nreads(nr))
              r = 'Inf';
            else
              r = num2str(Nreads(nr));
            end
            basicName = ['ss_','bac_dist_flag','_',num2str(bac_dist_flag),...
                         '_Nbac_in_mixture','_',num2str(Nbac_in_mixture(nb)),...
                         '_npower','_',num2str(npower*10),...
                         '_readlen','_',num2str(readlen),...
                         '_numIter','_',num2str(ii),...
                         '_Nreads','_',r];
            
            
            
            if npower==0.5
              np = '05';
              p = 1;
            else
              np = '0';
              p = 2;
            end
            load([userDir,'/CS/BAC/dataForSim/set/sec/',basicName,'_Nbacmix_',num2str(Nbac_in_mixture(nb)),'_Nread_',r,'_Readlen_',num2str(readlen),'_Npower_',np,'_bacdistflag_',num2str(bac_dist_flag),'_NoReads'],'correctWeight')
              
            
            try
              
              load([userDir,'/CS/BAC/sec/',basicName]);
              res(bac_dist_flag+1,nb,p,nr,ii) = max(abs(found{end})-correctWeight');
            catch
              basicName
            end
            
            
            
            %pause
            
         end
        
        end
      
      end
    end
  end
end


for bac_dist_flag=[0,1]
  for nb = 1:length(Nbac_in_mixture)
    for npower=[0.5 0]
      for readlen=[50]
        for nr = 1:length(Nreads)
          if npower==0.5
              np = '05';
              p = 1;
            else
              np = '0';
              p = 2;
          end
          r(bac_dist_flag+1,nb,p,nr) = max(res(bac_dist_flag+1,nb,p,nr,:));
        end
      end
    end
  end
end


% random bacteria;power law
disp('random bacteria;power law')
squeeze(r(1,:,1,:))
% random bacteria;uniform
disp('random bacteria;uniform')
squeeze(r(1,:,2,:))
% close bacteria; power law
disp('close bacteria; power law')
squeeze(r(2,:,1,:))
% close bacteria; uniform
disp('close bacteria; uniform')
squeeze(r(2,:,2,:))


