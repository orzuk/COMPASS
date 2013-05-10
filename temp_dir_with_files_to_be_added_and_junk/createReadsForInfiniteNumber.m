function [uniqueReads,uniqueReads_length,fracRelevantReadsForInfinity]=createReadsForInfiniteNumber(ind_bac_in_mix,correctWeight,readLength,basicSeqNameDir,basicSeqKey)

deleteDependentFlag = 0; % do not delete the dependent
[normalizedBac values independentInd,dependentInd] = prepareGroupOf1000DistributedSequenceFiles4(readLength,ind_bac_in_mix,basicSeqNameDir,basicSeqKey,deleteDependentFlag);

uniqueReads = values;
%keyboard

%keyboard
A = spalloc(size(normalizedBac,1),size(normalizedBac,2),nnz(normalizedBac));
A(find(normalizedBac)) =1;
uniqueReads_length = sum(A,2);

x = correctWeight(ind_bac_in_mix);
fracRelevantReadsForInfinity = normalizedBac*x;