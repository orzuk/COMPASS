function [normalizedBac values,varargout]=prepareMatrix(readLength,tmpInd,basicSeqNameDir,basicSeqKey)


% Sequence1 are the sequences of the current group
% load Sequence1 based on tmpInd
load(basicSeqKey,'positionInPart','len_uni')

numBAC = length(tmpInd);

Sequence1 = cell(numBAC,1);
%keyboard
for i=1:numBAC
  clear seq_*
  load([basicSeqNameDir,'seq_part_',num2str(positionInPart(tmpInd(i)))],['seq_',num2str(tmpInd(i))]);
  w = ['Sequence1{i} = seq_',num2str(tmpInd(i)),'{1};'];
  eval(w);
    
end
clear seq_* positionInPart

% seq is the list of possible sequences for each Sequence1
[normalizedBac values] = BuildMixingMatrixFromSequences(readLength,Sequence1,len_uni(tmpInd));

% check the rank
varargout{1} = Sequence1;
varargout{2} = len_uni(tmpInd);
