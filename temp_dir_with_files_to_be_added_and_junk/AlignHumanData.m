basicSeqNameDir = ['packed64/'];
basicSeqKey= ['keyNoNonACGT'];

disp('loading basicseqkey');
load(basicSeqKey,'len_uni');

NUMDBSEQS=410849;
USEREVERSECOMPLEMENT=0;
MAXIDEN=1000;

cppos=zeros(length(typseq),1);
ppos=zeros(MAXIDEN,length(typseq));
numiden=zeros(length(typseq),1);
disp('starting');
tic
for cseq=1:NUMDBSEQS
    if (mod(cseq,50)==0)
        toc
        disp(cseq);
        save('AlignHumanData-output','ppos','cseq');
        tic
    end
    [Header1,Sequence1] = loadSeqNames(cseq,basicSeqNameDir,basicSeqKey);
    for cshortseq=1:length(typseq)
        primerseq=nt2int(typseq{cshortseq});
        seq = unpack_seqs(Sequence1{1},len_uni(cseq) , 64);
        if (USEREVERSECOMPLEMENT==1)
            seq = seqrcomplement(seq);
        end
        [startpos,matchlen]=esubvector(primerseq,seq,200,400);
        if (matchlen>length(primerseq)-1)
            numiden(cshortseq)=numiden(cshortseq)+1;
            cppos(cshortseq)=cppos(cshortseq)+1;
            if (cppos(cshortseq)>MAXIDEN)
                cppos(cshortseq)=MAXIDEN;
            end
            ppos(cppos(cshortseq),cshortseq)=cseq;
        end
    end
end
