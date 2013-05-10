function prepareData_postProcessing(sampleName)

%keyboard
load(['~/CS/BAC/resSolexa_test_opt3_th50_mar25_2012/sol_test_opt3_th50_mar25_2012_',sampleName,'_noCorrection_100'])
load(['~/CS/BAC/12Samples/Solexa/data/opt3_th50_mar25_2012/readsSolexa_with_equalization_opt3_th50_mar25_2012_100_',sampleName])
load(['~/CS/BAC/12Samples/Solexa/data/opt3_th50_mar25_2012/illumina_reads100_',sampleName,'_unifreq_corr_opt3_th50'],'freq_read')

uniqueReads_length_no_equalization = freq_read;

readLength = 100;

userDir = getuserdir;
datasetName = 'primers750_primer_tail';
basicSeqNameDir = [userDir,'/CS/BAC/',datasetName,'/datNoNonACGT/packed64/'];
basicSeqKey = [userDir,'/CS/BAC/',datasetName,'/datNoNonACGT/keyNoNonACGT_',datasetName];

tmpInd = find(abs(found{end})>10^-3);
[normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOrFourth(readLength,tmpInd,basicSeqNameDir,basicSeqKey);
 
dataIn = struct;
[fracRelevantReads,sumRelevantReads] = currReadsFourth(uniqueReads,uniqueReads_length,values,0,dataIn);

[fracRelevantReads_no_equalization,sumRelevantReads] = currReadsFourth(uniqueReads,uniqueReads_length_no_equalization,values,0,dataIn);


[x_reg] = runOneGroupOf1000ForCompilationFourth(normalizedBac,fracRelevantReads);    

% check l1
[x_l1] = testL1(normalizedBac,fracRelevantReads);    

[x_l1_2_nor] = testL1_2(normalizedBac./max(fracRelevantReads),fracRelevantReads./max(fracRelevantReads));    
[x_l1_2_nonor] = testL1_2(normalizedBac,fracRelevantReads);    


N_mat = (fracRelevantReads*ones(1,size(normalizedBac,2))).*normalizedBac;
N_y = fracRelevantReads.^2;
mx = max(N_y);

N_y = N_y/mx;
N_mat = N_mat./mx;
[x_weighted] = runOneGroupOf1000ForCompilationFourth(N_mat,N_y);    

%keyboard
sqrt(sumsqr(normalizedBac*x_reg-fracRelevantReads))
sqrt(sumsqr(normalizedBac*x_l1-fracRelevantReads))
sqrt(sumsqr(N_mat*x_weighted-N_y))

subplot(1,5,1)
plot(x_reg,'.')
subplot(1,5,2)
plot(x_l1_2_nor,'.')
subplot(1,5,3)
title('l1 second')
plot(x_l1_2_nonor,'.')
title('l1 second')
subplot(1,5,4)
plot(x_l1,'.')
subplot(1,5,5)
plot(x_weighted,'.')
title(sampleName)




save(['~/CS/BAC/res12Samples_postprocessing/data_postProcessing_',sampleName],'fracRelevantReads','normalizedBac','found','tmpInd','values','uniqueReads','uniqueReads_length','uniqueReads_length_no_equalization','x_reg','x_l1','x_l1_2_nor','x_l1_2_nonor','x_weighted','fracRelevantReads_no_equalization','sumRelevantReads')        

%keyboard