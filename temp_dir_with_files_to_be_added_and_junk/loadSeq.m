function [Header1,Sequence1]=loadSeq(tmpInd,basicSeqNameDir,basicSeqKey)

load(basicSeqKey,'positionInPart','len_uni')

numBAC = length(tmpInd);
Header1 = cell(numBAC,1);
for i=1:numBAC
    
  load([basicSeqNameDir,'seq_part_',num2str(positionInPart(tmpInd(i)))],['seq_',num2str(tmpInd(i))]);
  w = ['Sequence1{i} = seq_',num2str(tmpInd(i)),'{1};'];
  eval(w);
  
  
end

