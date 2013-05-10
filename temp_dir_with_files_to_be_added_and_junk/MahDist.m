function [idenFracDist,idenFracDistMax,scoreDist]=MahDist(origSeqs,recSeqs,recFreq,threshold,auxData)
%keyboard

%keyboard
disp('Mahalanobis Distance');
disp(['loading ' num2str(length(origSeqs)) ' orig sequences']);

load(auxData.basicSeqKey,'len_uni')

%load the sequences
tmpInd = origSeqs;
[HeaderOrig,SequenceOrig] = loadSeqNames(tmpInd,auxData.basicSeqNameDir,auxData.basicSeqKey);
origCSeq=cell(1,length(tmpInd));
for a=1:length(tmpInd)
    origCSeq{a}=int2nt(unpack_seqs(SequenceOrig{a},len_uni(tmpInd(a)) , 64));
end

recFreq=abs(recFreq);
disp(['loading ' num2str(sum(recFreq>threshold)) ' recreated sequences']);
recSeqs=recSeqs(recFreq>threshold);
recFreq=recFreq(recFreq>threshold);
tmpInd = recSeqs;
[HeaderRec,SequenceRec] = loadSeqNames(tmpInd,auxData.basicSeqNameDir,auxData.basicSeqKey);
recCSeq=cell(1,length(tmpInd));
for a=1:length(tmpInd)
    recCSeq{a}=int2nt(unpack_seqs(SequenceRec{a},len_uni(tmpInd(a)) , 64));
end

disp('started comparison');

idenFracDist=zeros(length(origSeqs),length(recSeqs));
idenFracDistMax=zeros(length(origSeqs),length(recSeqs));
scoreDist=zeros(length(origSeqs),length(recSeqs));
for os=1:length(origSeqs)
    disp(os);
    seq1=origCSeq{os};
    tic
    for rs=1:length(recSeqs)
        seq2=recCSeq{rs};
        [Score, calgn, Start] = swalign(seq1, seq2, 'Alphabet','nt');
        idenFracDist(os,rs)= 1 - (sum(calgn(2,:)=='|') / min(length(seq1),length(seq2)));
        idenFracDistMax(os,rs)= 1 - (sum(calgn(2,:)=='|') / max(length(seq1),length(seq2)));
        scoreDist(os,rs)=Score;
    end
    toc
end
