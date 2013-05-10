
% create reads
clear

dirName = 'opt3_th50_mar25_2012';
typeName = '_unifreq_corr_opt3_th50';
datasetName = 'primers750_primer_tail';

readLength = 100;
createUniqueReads_12SamplesSolexa_with_eq_unicov_corrfreq(dirName,readLength,'O7',typeName)
createUniqueReads_12SamplesSolexa_with_eq_unicov_corrfreq(dirName,readLength,'O10',typeName)
createUniqueReads_12SamplesSolexa_with_eq_unicov_corrfreq(dirName,readLength,'S7',typeName)
createUniqueReads_12SamplesSolexa_with_eq_unicov_corrfreq(dirName,readLength,'S10',typeName)
createUniqueReads_12SamplesSolexa_with_eq_unicov_corrfreq(dirName,readLength,'M1',typeName)
createUniqueReads_12SamplesSolexa_with_eq_unicov_corrfreq(dirName,readLength,'M2',typeName)
createUniqueReads_12SamplesSolexa_with_eq_unicov_corrfreq(dirName,readLength,'M3',typeName)
createUniqueReads_12SamplesSolexa_with_eq_unicov_corrfreq(dirName,readLength,'M4',typeName)
createUniqueReads_12SamplesSolexa_with_eq_unicov_corrfreq(dirName,readLength,'S1',typeName)
createUniqueReads_12SamplesSolexa_with_eq_unicov_corrfreq(dirName,readLength,'S2',typeName)
createUniqueReads_12SamplesSolexa_with_eq_unicov_corrfreq(dirName,readLength,'S3',typeName)
createUniqueReads_12SamplesSolexa_with_eq_unicov_corrfreq(dirName,readLength,'S4',typeName)

% prepare the run files equalization

testData12Samples_with_eq_unicov_corrfreq_again({'S1','O7','O10','S7','S10','M1','M2','M3','M4','S2','S3','S4'},['test_',dirName],dirName,datasetName)
createFilesFor_testData12Samples_with_eq_unicov(['test_',dirName])



unix(['cd /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/;mkdir resSolexa_test_',dirName,';','rm -f /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/list_test_',dirName,';find sol_test_',dirName,'*  -name "*.mat" > list_',dirName])
% bash
w = ['for i in `cat list_',dirName,'`; do cp $i resSolexa_test_',dirName,';done'];


% create the matrix
clear
load ~/CS/BAC/resSolexa_test_opt3_th50_mar25_2012/sol_test_opt3_th50_mar25_2012_S1_noCorrection_100
load /home/csfaculty/shental/CS/BAC/12Samples/Solexa/data/opt3_th50_mar25_2012/readsSolexa_with_equalization_opt3_th50_mar25_2012_100_S1

readLength = 100;

userDir = getuserdir;
datasetName = 'primers750_primer_tail';
basicSeqNameDir = [userDir,'/CS/BAC/',datasetName,'/datNoNonACGT/packed64/'];
basicSeqKey = [userDir,'/CS/BAC/',datasetName,'/datNoNonACGT/keyNoNonACGT_',datasetName];

tmpInd = find(abs(found{3})>10^-3);
[normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOrFourth(readLength,tmpInd,basicSeqNameDir,basicSeqKey);
  
  
dataIn = struct;
[fracRelevantReads,sumRelevantReads] = currReadsFourth(uniqueReads,uniqueReads_length,values,0,dataIn);

[x0] = runOneGroupOf1000ForCompilationFourth(normalizedBac,fracRelevantReads);    

clear tk iglist*
load ~/iglistall
tk = 1:size(normalizedBac,1);
tk(iglistall) = [];

[x1] = runOneGroupOf1000ForCompilationFourth(normalizedBac(tk,:),fracRelevantReads(tk));    
save ~/solution_S1_without_iglistall x

clear tk iglist*
load ~/iglistall-alot
tk = 1:size(normalizedBac,1);
tk(iglistall) = [];

[x2] = runOneGroupOf1000ForCompilationFourth(normalizedBac(tk,:),fracRelevantReads(tk));    
save ~/solution_S1_without_iglistall_alot x

clear tk iglist*
load ~/iglistall-wobble
tk = 1:size(normalizedBac,1);
tk(iglistall) = [];

[x3] = runOneGroupOf1000ForCompilationFourth(normalizedBac(tk,:),fracRelevantReads(tk));    
save ~/solution_S1_without_iglistall_wobble x

load('~/CS/BAC/12Samples/Solexa/data/opt3_th50_mar25_2012/illumina_reads100_S1_unifreq_corr_opt3_th50','freq_read')
new_uniqueReads_length = freq_read;
[new_fracRelevantReads,new_sumRelevantReads] = currReadsFourth(uniqueReads,new_uniqueReads_length,values,0,dataIn);

save ~/new_fracRelevantReads new_fracRelevantReads




a = sqrt(sumsqr((normalizedBac*x-fracRelevantReads)))

b = sqrt(sumsqr((normalizedBac(:,[1 ])*x([1])-fracRelevantReads)))


% check l1
[x] = testL1(normalizedBac,fracRelevantReads);    
%x = x./sum(x);

[x1_2] = testL1_2(normalizedBac./max(fracRelevantReads),fracRelevantReads./max(fracRelevantReads));    
sqrt(sumsqr(normalizedBac*x1_2-fracRelevantReads))


N_mat = (fracRelevantReads*ones(1,size(normalizedBac,2))).*normalizedBac;
N_y = fracRelevantReads.^2;
mx = max(N_y);

N_y = N_y/mx;
N_mat = N_mat./mx;
[x_weighted] = runOneGroupOf1000ForCompilationFourth(N_mat,N_y);    

save ~/CS/BAC/x_weighted x_weighted


save ~/CS/BAC/dataForS1_opt3 fracRelevantReads found tmpInd normalizedBac values uniqueReads uniqueReads_length


%%%%%%%%%%%%%%%%%555
% data preprocessing

samples = {'S1','O7','O10','S7','S10','M1','M2','M3','M4','S2','S3','S4'};

for i=1:length(samples)
  figure(i)
  prepareData_postProcessing(samples{i});
end

