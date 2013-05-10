function [uniqueReads,uniqueReads_length,fracRelevantReadsForInfinity]=createReadsForInfiniteNumberOr(ind_bac_in_mix,correctWeight,readLength,basicSeqNameDir,basicSeqKey)

[normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOr(readLength,ind_bac_in_mix,basicSeqNameDir,basicSeqKey);

uniqueReads = values;
%keyboard

%keyboard
A = spalloc(size(normalizedBac,1),size(normalizedBac,2),nnz(normalizedBac));
A(find(normalizedBac)) =1;
uniqueReads_length = sum(A,2);

x = correctWeight(ind_bac_in_mix);
fracRelevantReadsForInfinity = normalizedBac*x;