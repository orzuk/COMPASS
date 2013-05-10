function collectThoseThatdidNotRun(basicPrefix,outFileDirName,numProcessors,bdf,nb,np,numIter,nr)
disp('works with non-ACGT data')
%userDir = getuserdir;
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
load([userDir,'/CS/BAC/cellsimbac2'])


fid = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','problamatic_',basicPrefix,'_iter_',num2str(iterationNumber)],'r');
list = textscan(fid,'%d %s');
fclose(fid);


createReadsFlag = 0;
N = 410849;
numRepeatRandomGroups = 10;

% find those that were not run
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
            
             
          end
        end
      end
    end
  end
end




unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','listOfR_',basicPrefix,'*']);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','run_*'])
save([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',basicPrefix,'_list'],'totalListOfNames','totalListOfRuns')
  

