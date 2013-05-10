function [res,resCell]=loadResOr(basicPrefix,outFileDirName,bdf,nb,np,numIter,nr)
disp('works with non-ACGT data')
%userDir = getuserdir;
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';

r = [];

createReadsFlag = 0;
N = 410849;
numRepeatRandomGroups = 10;

for bac_dist_flag=bdf
  for nbb = 1:length(nb)
    for npower=np
      for readlen=[50]
        for nrr=1:length(nr)
          for ii=1:numIter
         
            
            Nbac_in_mixture = nb(nbb);
            Nreads = nr(nrr);
            
            if isinf(Nreads)
              rr = 'Inf';
            else
              rr = num2str(Nreads);
            end
            
            if npower==0.5
              curr_np = '05';
              p = 1;
            else
              curr_np = '0';
              p = 2;
            end
            
            basicName = [basicPrefix,'_bac_dist_flag','_',num2str(bac_dist_flag),...
                         '_Nbac_in_mixture','_',num2str(Nbac_in_mixture),...
                         '_npower','_',num2str(npower*10),...
                         '_readlen','_',num2str(readlen),...
                         '_numIter','_',num2str(ii),...
                         '_Nreads','_',num2str(Nreads)];
            
            load([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',basicName,'_Nbacmix_',num2str(Nbac_in_mixture),'_Nread_',rr,'_Readlen_',num2str(readlen),'_Npower_',curr_np,'_bacdistflag_',num2str(bac_dist_flag),'_NoReads'],'correctWeight');
            
            try
              
              load([userDir,'/CS/BAC/',outFileDirName,'/',basicName]);
              
              %keyboard
              nonEmpty = find(cellfun(@isempty,found)==0);
              tmpRes = sparse(1+length(nonEmpty),length(correctWeight));
              nonZer = find(correctWeight);
              
              tmpRes = cell(1+length(nonEmpty),1);
              tmpRes{1} = [nonZer,correctWeight(nonZer)];
              for zz=1:length(nonEmpty)
                nonZer = find(abs(found{nonEmpty(zz)})>10^-5);
                tmpRes{zz+1} = [nonZer',found{nonEmpty(zz)}(nonZer)'];
              end
              
              
              
              %tmpRes(1,nonZer) = correctWeight(nonZer);
              
              %for zz=1:length(nonEmpty)
              %  nonZer = find(abs(found{nonEmpty(zz)})>10^-4);
              %  tmpRes(zz+1,nonZer) = found{nonEmpty(zz)}(nonZer);
              %end
              resCell{bac_dist_flag+1,nbb,p,nrr,ii} = tmpRes;
              % resCell(bac_dist_flag+1,nbb,p,nrr,ii) = tmpRes;
              
              res{bac_dist_flag+1,nbb,p,nrr,ii} = max(abs(found{end})-correctWeight');
            catch
              basicName
            end
            
            
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
            
            
            
            
            
            
          end % num iter
          %r(bac_dist_flag+1,nbb,p,nrr) = max(res(bac_dist_flag+1,nbb,p,nrr,:));
        end
      end
    end
  end
end
