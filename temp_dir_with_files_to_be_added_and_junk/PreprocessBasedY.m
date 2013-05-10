% input :
% mat_i_j = the fraction of times k-mer i appears in species j
% reads_i - the fraction of reads which correspond to k-mer i
% output:
% fracReads_i - the fraction of k-mers of species i which are present in
% the reads
function [fracReads]=PreprocessBasedY(mat,reads)
numOfSeqs=size(mat,2);
fracReads=zeros(numOfSeqs,1);
for a=1:numOfSeqs
    merIsPresent=find(mat(:,a)>0);
    cdat=reads(merIsPresent);
    fracReads(a)=length(find(cdat>0))/length(merIsPresent); % fraction of totla k-mers present
end
