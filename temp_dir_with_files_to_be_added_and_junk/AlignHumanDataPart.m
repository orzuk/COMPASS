function AlignHumanDataPart(i,numiden,cppos,ppos,part,typseq,len_uni,USEREVERSECOMPLEMENT,auxData)

MAXIDEN = size(ppos,1);

for cseq=part(1):part(2)-1

   if (mod(cseq,50)==0)
     tic
       disp(cseq);
       save(['~/CS/BAC/AlignHumanData_output_',num2str(i)],'ppos','cseq');
     toc
   end
   [Header1,Sequence1] = loadSeqNames(cseq,auxData.basicSeqNameDir,auxData.basicSeqKey);
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
