function [dist,rmse,maxdist,numrecreated,maxnonfreq,ingroupmaxdist,rec,prec,thresh]=compareResultsCellOr(origInd,origWeight,recInd,recWeight,basicSeqKey,basicSeqNameDir)

if 1==2
basicSeqKey = '~/CS/BAC/datNoNonACGT/keyNoNonACGT';%'E:\Amnon\NGSeqCS\CompareSolution\dat450000\key450000';
basicSeqNameDir = '~/CS/BAC/datNoNonACGT/packed64/';%'E:\Amnon\NGSeqCS\CompareSolution\dat450000\';
numbact=4;
distrib=2;
nreads=3;
cf=1;

load ~/CS/BAC/res_t4_bac_dist_flag_1_0_Nbac_in_mixture_10_50_100_500_npower_0_5_readlen_10000_100000_1000000_Inf_numIter_10_Nreads_10000_100000_1000000_Inf.mat

crun=resCell{1,numbact,distrib,nreads,cf};
            orig=crun{1};
            origInd=orig(:,1);
            origWeight=abs(orig(:,2));
            rec=crun{2}; %NEED to add different run stop condition stats
            recInd=rec(:,1);
            recWeight=abs(rec(:,2));
            [dist,rmse,maxdist,numrecreated,maxnonfreq,ingroupmaxdist,rec,prec,thresh]=compareResultsCellOr(origInd,origWeight,recInd,recWeight,basicSeqKey,basicSeqNameDir);
end


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
%keyboard
% Noam changes
[normalizedBac values Sequence1 len_uni1] = prepareGroupOf1000DistributedSequenceFilesOr(readLength,b,basicSeqNameDir,basicSeqKey); 

unpackedSequence1 = cell(length(Sequence1),1);
for i=1:length(Sequence1)
  unpackedSequence1{i} = int2nt(unpack_seqs(Sequence1{i},len_uni1(i),64)); % if you want to unpack it
end
% unpack Sequence1

%%%%% end Noam changes



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
    [freqset,igmd]=Compare2(normalizedBac,i_c1,origWeight,i_f1,recWeight,a,unpackedSequence1);
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