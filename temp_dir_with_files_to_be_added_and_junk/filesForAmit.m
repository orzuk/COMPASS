
1) load reads file

2) convert your reads to packed format

uniqueReads = pack_seqs(uni_read_seq,64);
uniqueReads_length = freq_corr2pos';

3) create the matrix and relevant vector

load('keyNoNonACGT_primers750_primer_tail','positionInPart','len_uni') % len_uni
readLength = 100;
tmpInd = are the indices of the bacteria

[normalizedBac values] = BuildMixingMatrixFromSequences(readLength,your_Sequence,len_uni(tmpInd));


[fracRelevantReads,sumRelevantReads] = currReadsFourth(uniqueReads,uniqueReads_length,values,0,struct);
4) solve
[x0] = runOneGroupOf1000ForCompilationFourth(normalizedBac,fracRelevantReads);    




