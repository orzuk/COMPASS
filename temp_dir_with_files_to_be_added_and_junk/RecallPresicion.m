% Calculate the weighted recall/precision of the reconstruction
% input:
% distmat - the distance matrix of all the original and reconstructed
% sequences (from mothur and squareform)
% origFreq - frequencies of original seqs (0 if not present)
% recFreq - frequencies of reconstructed seqs (0 if not present)
% simThresh - the sequence similarity threshold to define sequences as
% identical for the reconstruction
% output:
% rec,prec - recall (fraction of original seqs) and precision
% (fraction of correct seqs) of the reconstruction.
% need the reconstructed frequency to be
% 0.5-2 times the original frequency or <0.2% freq diff
% freq - the L2 distance of the frequency vectors
function [rec,prec,freq]=RecallPresicion(distmat,origFreq,recFreq,simThresh)
origFreq(origFreq==0)=-1;
recFreq(recFreq==0)=-1;
    recPresent=find(recFreq>=0);
    origPresent=find(origFreq>=0);
    
    totrec=zeros(length(origPresent),1);
    totfreq=zeros(length(origPresent),1);
    isrecok=zeros(length(recPresent),1);
    % variables for frequency threshold
    isrecokfreq=zeros(length(recPresent),1);
    isorigokfreq=zeros(length(origPresent),1);
    FTHRESH=1.2;
    FATHRESH=0.002;
    
    freq=0;
    
    [neardist,nearpos]=min(distmat(origPresent,recPresent));
    % go over all original sequences
    for cseq=1:length(origPresent)
        corig=origPresent(cseq);
        % find corresponding reconstructed seqs
        corseqs=find(nearpos==cseq);
        % select only ones similar enough
        corseqs=corseqs(neardist(corseqs)<=simThresh);
        if (~isempty(corseqs))
            isrecok(corseqs)=1;
            ofreq=origFreq(corig);
            rfreq=sum(recFreq(recPresent(corseqs)));
            if (((ofreq/rfreq<FTHRESH)&&(ofreq/rfreq>1/FTHRESH))||(abs(ofreq-rfreq)<FATHRESH))
                freq=freq+(ofreq-rfreq)^2;
%                freq=freq+abs(ofreq-rfreq)/ofreq;
                isrecokfreq(corseqs)=1;
                isorigokfreq(cseq)=1;
            else
                freq=freq+ofreq^2;
            end
        end
        totrec(cseq)=length(corseqs);
        totfreq(cseq)=sum(recFreq(recPresent(corseqs)));
    end
%    trec=sum(totrec>0)/length(origPresent);
%    tprec=sum(totrec>0)/length(recPresent);
%    tprec=sum(recFreq(recPresent(isrecok>0)))/sum(recFreq(recPresent));
%    trec=sum(origFreq(origPresent(totrec>0)));

%    tprec=sum(recFreq(recPresent(isrecok>0)))/sum(recFreq(recPresent));
%    trec=sum(origFreq(origPresent(totrec>0)));
    tprec=sum(recFreq(recPresent(isrecokfreq>0)))/sum(recFreq(recPresent));
    trec=sum(origFreq(origPresent(isorigokfreq>0)));
    freq=freq+sum(recFreq(recPresent(isrecokfreq==0)).^2);
%    freq=freq/length(origPresent(isorigokfreq>0));
    rec=trec;
    prec=tprec;