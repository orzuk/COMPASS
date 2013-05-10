function res=extractRes(dirName,basicName,bdf,nb,np,nr,addNoise,correctMeasurementForReadErrors,readLength,N,howFlag)
%keyboard
%keyboard
dr = dir([dirName,basicName,'_bac_dist_flag_',num2str(bdf),'_Nbac_in_mixture_',num2str(nb),'_npower_',num2str(np*10),'_readlen_',num2str(readLength),'_*','_Nreads_',num2str(nr),'_Noise_',num2str(addNoise),'_Correction_',num2str(correctMeasurementForReadErrors),'.mat']);
tit = [dirName,basicName,'_bac_dist_flag_',num2str(bdf),'_Nbac_in_mixture_',num2str(nb),'_npower_',num2str(np*10),'_readlen_',num2str(readLength),'_*','_Nreads_',num2str(nr),'_Noise_',num2str(addNoise),'_Correction_',num2str(correctMeasurementForReadErrors),'.mat'];

if ~exist('howFlag','var')
  howFlag = 0;
end
  
res = NaN*ones(1,length(dr));
if howFlag
  for i=1:length(dr)
    r = load([dirName,dr(i).name]);
    
    x_found = zeros(1,N);
    x_found(r.resCell{2}(:,1)) = abs(r.resCell{2}(:,2));
    
    x_correct = zeros(1,N);
    x_correct(r.resCell{1}(:,1)) = r.resCell{1}(:,2);
  
    res(i) = max(abs(x_found-x_correct));
    
  end
else
  for i=1:length(dr)
    r = load([dirName,dr(i).name]);
    x_found = abs(r.found{length(r.found)});
    res(i) = max(abs(x_found-r.correctWeight'));
    
  end
end
%keyboard

disp('should consider not 2 extractRes.m !!!')
max(res)
[readLength length(find(res>1/nb-10^-4))]
figure;
plot(res,'.')
title(['bdf: ',num2str(bdf),' nb: ',num2str(nb),' np: ',num2str(np*10),' readlen: ',num2str(readLength),' Nreads:',num2str(nr),' Noise: ',num2str(addNoise),' Correction:',num2str(correctMeasurementForReadErrors)])


%keyboard
res = [mean(res), std(res)];


%plot(x_found,x_correct,'.')


if 1==2
  for j=1:length(r.save_store_kp)
    setdiff(r.resCell{1}(:,1),r.save_store_kp{j})
  end

  basicSeqNameDir = ['~/CS/BAC/datNoNonACGT/packed64/'];
  basicSeqKey= ['~/CS/BAC/datNoNonACGT/keyNoNonACGT'];
  
  load(basicSeqKey,'len_uni') 

  z = find(x_found>0.003 & x_correct==0)
  z_c = find(x_found<0.006 & x_correct==0.01)
  
  tmpZ = z;
  z = z(1);
  z_c = tmpZ(2);
%  z_c = z_c(2);
  
  
  [Header1,Sequence1] = loadSeqNames(z,basicSeqNameDir,basicSeqKey);
  seq = int2nt(unpack_seqs(Sequence1{1},len_uni(z),64));
  
  
  
  clear ham
  for j=1:length(z_c)
    [Header1,Sequence_c] = loadSeqNames(z_c(j),basicSeqNameDir,basicSeqKey);
    seq_c = int2nt(unpack_seqs(Sequence_c{1},len_uni(z_c(j)),64));
    [cscore,algn]=swalign(seq,seq_c,'Alphabet','NT');
    ham(j) = length(find(algn=='|'));
  end
  
  ind_bac_in_mix = r.resCell{1}(:,1)';
  correctWeight = zeros(410849,1);correctWeight(ind_bac_in_mix) = r.resCell{1}(:,2);
  % output the raw data
  readLength = 50;
  [uniqueReads,uniqueReads_length,fracRelevantReadsForInfinity,normalizedBac_junk,values_junk] = createReadsForInfiniteNumberOrFourth(ind_bac_in_mix,correctWeight,readLength,basicSeqNameDir,basicSeqKey);
  % follows the order in BuildMixingMatrixFromSequences 
  
  %tmpInd = dr.resCell{2}(:,1)';
  tmpInd = ind_bac_in_mix;
  [normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOrFourth(readLength,tmpInd,basicSeqNameDir,basicSeqKey);
  % follows the order in BuildMixingMatrixFromSequences 
  
  [normalizedBac1 values1] = prepareGroupOf1000DistributedSequenceFilesOrFourth(readLength,tmpInd,basicSeqNameDir,basicSeqKey);
  
  %[normalizedBac1 values1] = prepareGroupOf1000DistributedSequenceFilesOrFourth(readLength,tmpInd,basicSeqNameDir,basicSeqKey);
  
  %prepareGroupOf1000DistributedSequenceFilesOrFourth(readLength,tmpInd,basicSeqNameDir,basicSeqKey);
  
  dataIn = struct;
  dataIn.fracRelevantReadsForInfinity = fracRelevantReadsForInfinity;
  [fracRelevantReads,junk] = currReadsFourth(uniqueReads,uniqueReads_length,values,1,dataIn);
  
  [x]=runOneGroupOf1000ForCompilationFourth(normalizedBac,fracRelevantReads);    
  
  new_found = zeros(1,410849);
  new_found(tmpInd) = x;
  
  
  
end
