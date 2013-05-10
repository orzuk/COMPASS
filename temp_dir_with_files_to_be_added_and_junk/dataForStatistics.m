basicSeqNameDir = ['~/CS/BAC/datNoNonACGT/packed64/'];
basicSeqKey= ['~/CS/BAC/datNoNonACGT/keyNoNonACGT'];

load(basicSeqKey,'len_uni') 

tmpInd = 1:1000;
[Header1,Sequence1] = loadSeqNames(tmpInd,basicSeqNameDir,basicSeqKey);

seq = cell(1,length(tmpInd));
for i=1:length(tmpInd)
  seq{i} = int2nt(unpack_seqs(Sequence1{i},len_uni(tmpInd(i)) , 64));
end


