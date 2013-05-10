function [dist,rmse,maxdist,numrecreated,maxnonfreq,ingroupmaxdist,rec,prec,thresh]=compareResultsCell(origInd,origWeight,recInd,recWeight,basicSeqKey,basicSeqNameDir)
  
%load(resFile)
%load(correctFile,'correctWeight','ind_bac_in_mix')

%disp(['using the last iteration which consists of ',num2str(length(store_kp{end})),' bacteria. If that is too small - use preceding iterations. see compareResults.m'])
%disp('thresholds the found frequency at 10^-3. see compareResults.m')

%solution = abs(found{end});
%foundBAC = find(abs(found{end})>10^-3)';

b=[origInd;recInd];
%b = [ind_bac_in_mix;foundBAC];
b = unique(b);
  
readLength = 50;
[normalizedBac values Sequence1] = prepareGroupOf1000DistributedSequenceFiles2(readLength,b,basicSeqNameDir,basicSeqKey); 

%[junk,i_f1] = ismember(foundBAC,b);
%[junk,i_c1] = ismember(ind_bac_in_mix,b);
[junk,i_f1] = ismember(recInd,b);
[junk,i_c1] = ismember(origInd,b);

dist=[];
rmse=[];
maxdist=[];
numrecreated=[];
maxnonfreq=[];
ingroupmaxdist=[];
cdist=1;
for a=1:100:1
    dist=[dist a];
    [freqset,igmd]=Compare2(normalizedBac,i_c1,origWeight,i_f1,recWeight,a,Sequence1);
    ingroupmaxdist=[ingroupmaxdist igmd];
    rmse=[rmse sqrt(sum((freqset(:,1)-freqset(:,2)).^2))];
    maxdist=[maxdist max(abs(freqset(:,1)-freqset(:,2)))];
    numrecreated=[numrecreated sum((freqset(:,1)>0).*(freqset(:,2)>0))];
    mnf=max([0;freqset(find(freqset(:,1)==0),2)]);
    maxnonfreq=[maxnonfreq mnf];
    for cthreshi=1:500
       cthresh=cthreshi*max(abs(freqset(:,1)))/500; 
        thresh(cthreshi)=cthresh;
        [trec,tprec]=RecallPrecision(freqset(:,1),freqset(:,2),cthresh);
        rec(cdist,cthreshi)=trec;
        prec(cdist,cthreshi)=tprec;
    end
    cdist=cdist+1;
%figure;PlotFreqs(freqset,dist);
end