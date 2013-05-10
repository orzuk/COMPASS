function r=loadRes(basicPrefix,bdf,nb,np,numIter,nr)
r = [];

userDir = getuserdir;
%userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
load([userDir,'/CS/BAC/cellsimbac2'])

createReadsFlag = 0;
N = 455055;
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
              np = '05';
              p = 1;
            else
              np = '0';
              p = 2;
            end
            
            
            basicName = [basicPrefix,'_bac_dist_flag','_',num2str(bac_dist_flag),...
                         '_Nbac_in_mixture','_',num2str(Nbac_in_mixture),...
                         '_npower','_',num2str(npower*10),...
                         '_readlen','_',num2str(readlen),...
                         '_numIter','_',num2str(ii),...
                         '_Nreads','_',num2str(Nreads)];
            
            
            load([userDir,'/CS/BAC/dataForSim/set/sec/',basicName,'_Nbacmix_',num2str(Nbac_in_mixture),'_Nread_',rr,'_Readlen_',num2str(readlen),'_Npower_',np,'_bacdistflag_',num2str(bac_dist_flag),'_NoReads'],'correctWeight');
            
            try
              
              load([userDir,'/CS/BAC/sec/',basicName]);
              
              res(bac_dist_flag+1,nbb,p,nrr,ii) = max(abs(found{end})-correctWeight');
            catch
              basicName
            end
            
          end % num iter
          r(bac_dist_flag+1,nbb,p,nrr) = max(res(bac_dist_flag+1,nbb,p,nrr,:));
          
        end
      end
    end
  end
end


res