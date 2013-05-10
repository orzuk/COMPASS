function [Header1,Sequence1]=loadSeqNames(tmpInd,basicSeqNameDir,basicSeqKey)

load(basicSeqKey,'positionInPart','len_uni')

numBAC = length(tmpInd);
Header1 = cell(numBAC,1);
for i=1:numBAC
  clear head_*
  load([basicSeqNameDir,'head_part_',num2str(positionInPart(tmpInd(i)))],['head_',num2str(tmpInd(i))]);
  w = ['Header1{i} = head_',num2str(tmpInd(i)),'{1};'];
  eval(w);
    
  load([basicSeqNameDir,'seq_part_',num2str(positionInPart(tmpInd(i)))],['seq_',num2str(tmpInd(i))]);
  w = ['Sequence1{i} = seq_',num2str(tmpInd(i)),'{1};'];
  eval(w);
  
  
end

