clear 
load ~/CS/BAC/canonicalFormWithL1InLastStage/tmp


bac = randperm(size(normalizedBac,2));
bac = bac(1:10:end);


[x1_2_nor] = testL1_2(normalizedBac(:,bac)./max(fracRelevantReads),fracRelevantReads./max(fracRelevantReads));    

opts = optimset('linprog');
set( opts, 'TolFun', 1e-12 );

[x1_2_m,resnorm,res,ex] = madlin(normalizedBac(:,bac)./max(fracRelevantReads),fracRelevantReads(:)./max(fracRelevantReads),[],[],[],[],zeros(1,length(bac)),[],[],opts);

x1_lp = linprog( normalizedBac./max(fracRelevantReads), [], [], fracRelevantReads./max(fracRelevantReads), y, lb, [], [], opts );

plot(x1_2_nor,x1_2_m,'.')


clear

%tmpInd = [1:10^3:10^5-1];
tmpInd = (1:500);
basicSeqNameDir = '~/CS/BAC/datNoNonACGT/packed64/';
basicSeqKey = '~/CS/BAC/datNoNonACGT/keyNoNonACGT';
[normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOrFourth(100,tmpInd,basicSeqNameDir,basicSeqKey);
  
%c = zeros(length(tmpInd),1);
%c = (1/length(c))*ones(length(tmpInd),1);
c = rand(1,length(tmpInd));
c = c./sum(c);
c = c';
plot(c,'.')


fracRelevantReads = normalizedBac*c;
fracRelevantReads1 = fracRelevantReads.*(1+1/2*abs(randn(length(fracRelevantReads),1)));

a = 1000;
fracRelevantReads1 = floor(fracRelevantReads*a)/a;

%a = randperm(length(fracRelevantReads));
%fracRelevantReads(a(1:0.6*length(fracRelevantReads))) = 0;
[x1_2_nor] = testL1_2(normalizedBac,fracRelevantReads1);    
[x2_2_nor] = testL2_2(normalizedBac,fracRelevantReads1);    

figure(1)
plot(c,x1_2_nor,'.')
figure(2)
plot(c,x2_2_nor,'.')





%%%%%%%%%%%%%55
clear
load ~/CS/BAC/canonicalFormWithL1InLastStage/canonicalFormWithL1InLastStage_bac_dist_flag_0_Nbac_in_mixture_1000_npower_5_readlen_100_numIter_1_Nreads_1000000_Noise_0_Correction_0.mat
load ~/CS/BAC/canonicalFormWithL1InLastStage/canonicalFormWithL1InLastStage_bac_dist_flag_0_Nbac_in_mixture_1000_npower_5_readlen_100_numIter_1_Nreads_1000000_Noise_0_Correction_0_L1data.mat

basicSeqNameDir = '~/CS/BAC/datNoNonACGT/packed64/';
basicSeqKey = '~/CS/BAC/datNoNonACGT/keyNoNonACGT';
[normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOrFourth(100,tmpInd,basicSeqNameDir,basicSeqKey);
  
dataIn = struct;
[fracRelevantReads,sumRelevantReads] = currReadsFourth(uniqueReads,uniqueReads_length,values,0,dataIn);
%fracRelevantReads = myCorrectReadsNew3_for_SL4(values,uniqueReads,uniqueReads_length');

[x1_2_nor] = testL1_2(normalizedBac,fracRelevantReads);    


plot(x1_2_nor,'.')

10e6 10e7
