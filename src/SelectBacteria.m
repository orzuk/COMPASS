% Returns a reduced set of bacteria (hopefully without removing true ones)
% Input
% bactnat - the bacteria matrix (#k-mers X #species)
% reads - the next gen read frequency for each k-mer
% percentile - the percentile cutoff to use (try 20)
%
% Output
% outset - the ids of the bacteria to keep
function [outset]=SelectBacteria(bactmat,reads,percentile)
numOfSeqs=size(bactmat,2);
minFreq=zeros(numOfSeqs,1);
for a=1:numOfSeqs
    merIsPresent=find(bactmat(:,a)>0);
    cdat=reads(merIsPresent);
    minFreq(a)=prctile(cdat,percentile);
end
outset=find(minFreq>0);
end