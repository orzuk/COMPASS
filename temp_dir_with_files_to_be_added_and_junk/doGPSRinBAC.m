function doGPSRinBAC(uniqueReads,uniqueReads_length,tmpInd,auxData)


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
M0 = normalizedBac;
s = sum(M0,1);
M0 = M0./(ones(size(M0,1),1)*s);

qMeasurement = fracRelevantReads./sum(fracRelevantReads);

tau = 0.0005*max(abs(M0'*qMeasurement));
fractionalOutput = applyGPSR(qMeasurement,M0,tau);
