function x1=doCSinBAC(uniqueReads,uniqueReads_length,tmpInd,auxData)


if auxData.brFlag==1 &  ~isempty(findstr(getenv('HOSTNAME'),'openu'))
  userDir = getuserdir;
  auxData.basicSeqNameDir = strrep(auxData.basicSeqNameDir,'/seq/orzuk2/compressed_sensing/metagenomics/next_gen',userDir);
  auxData.basicSeqKey = strrep(auxData.basicSeqKey,'/seq/orzuk2/compressed_sensing/metagenomics/next_gen',userDir);
end

[normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOrFourth(auxData.readLength,tmpInd,auxData.basicSeqNameDir,auxData.basicSeqKey);
  
dataIn = struct;
[fracRelevantReads,sumRelevantReads] = currReadsFourth(uniqueReads,uniqueReads_length,values,0,dataIn);

if auxData.addNoiseFlag==1

  if auxData.readLength>=32 & auxData.readLength<64
    fracRelevantReads = myCorrectReadsNew3_for_SL2(values,uniqueReads,uniqueReads_length');
  elseif auxData.readLength>=64 & auxData.readLength<128
    fracRelevantReads = myCorrectReadsNew3_for_SL4(values,uniqueReads,uniqueReads_length');
    
  elseif auxData.readLength>=128 & auxData.readLength<256
    fracRelevantReads = myCorrectReadsNew3_for_SL8(values,uniqueReads,uniqueReads_length');
  elseif auxData.readLength>=256 & auxData.readLength<384
    fracRelevantReads = myCorrectReadsNew3_for_SL12(values,uniqueReads,uniqueReads_length');
  elseif auxData.readLength>=384 & auxData.readLength<512
    fracRelevantReads = myCorrectReadsNew3_for_SL13(values,uniqueReads,uniqueReads_length');
  else 
    disp('should prepare myCorrectReadsNew for this read length. solveForGroupOrFourth.m')
  end

end

keyboard

k = randperm(size(normalizedBac,2));
k = k(1:2000);
k = setdiff(k,i2);
x1 = testL1_CS(normalizedBac(:,[i2,k]),fracRelevantReads,0);

M0 = normalizedBac;
s = sum(M0,1);
M0 = M0./(ones(size(M0,1),1)*s);
qMeasurement = fracRelevantReads./sum(fracRelevantReads);

[x_L1,x_L2,x_L1_and_L1x,x_L2_and_L1x]=testVarious(M0,qMeasurement,1);


x1 = testL1_CS(M0,qMeasurement);

run_lin_prog = 1;
A = full(M0'*M0);
[unique_inds x_min x_max] = set_solution_bounds(A, x1, qMeasurement, 0, run_lin_prog,10^-3);
      


s2 = s(find(x1>10^-3));
x1_a = testL1_CS(M0(:,s2),qMeasurement);

plot(s,x1,'.');hold on;
plot(s2,x1_a,'ro')

k = randperm(size(normalizedBac,2));
k = k(1:500);
k = setdiff(k,i2);

[x_L1,x_L2,x_L1_and_L1x,x_L2_and_L1x]=testVarious(M0(:,[i2(1:50),k]),qMeasurement,1);

subplot(1,4,1);
plot(x_L1,'.');
x = x_L1;title([norm(M0(:,[i2(1:50),k])*x-qMeasurement,1),norm(x,1)]);
subplot(1,4,2);
plot(x_L2,'.');
x = x_L2;title([norm(M0(:,[i2(1:50),k])*x-qMeasurement,1),norm(x,1)]);

subplot(1,4,3);
plot(x_L1_and_L1x,'.');
x = x_L1_and_L1x;title([norm(M0(:,[i2(1:50),k])*x-qMeasurement,1),norm(x,1)]);

subplot(1,4,4);
plot(x_L2_and_L1x,'.');
x = x_L2_and_L1x;title([norm(M0(:,[i2(1:50),k])*x-qMeasurement,1),norm(x,1)]);

%%%%%%%%%%%%%%%%%%%%%%%%%

tmpInd = unique([curr_kp{1},ind_bac_in_mix]);
[normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOrFourth(auxData.readLength,tmpInd,auxData.basicSeqNameDir,auxData.basicSeqKey);
  dataIn = struct;
[fracRelevantReads,sumRelevantReads] = currReadsFourth(uniqueReads,uniqueReads_length,values,0,dataIn);


M0 = normalizedBac;
%s = sum(M0,1);
%M0 = M0./(ones(size(M0,1),1)*s);
qMeasurement = fracRelevantReads;

[x_L1,x_L2,x_L1_and_L1x,x_L2_and_L1x]=testVarious(M0,qMeasurement,1);

sqrt(sumsqr(M0*x_L1-qMeasurement))
sqrt(sumsqr(M0*x_L2-qMeasurement))

sum(abs(M0*x_L1-qMeasurement))
sum(abs(M0*x_L2-qMeasurement))


c = zeros(size(M0,2),1);
[junk,i1,i2] = intersect(ind_bac_in_mix,tmpInd);
c(i2) = 1;

sqrt(sumsqr(M0*c-qMeasurement))


%%%%%%%%%%%%%%%%%


norm(M0*x-qMeasurement,1),norm(x,1)