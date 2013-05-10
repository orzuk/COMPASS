function [dat_f,dat_r]=compareMichael_forward_and_reverse(sampleNum,dataBaseFile)
%keyboard
bcs_data_path = '~/CS/BAC/';

load(dataBaseFile); 

clear fracRelevantReads sumRelevantReads reads F tmpMat
reads = cell(6,1);
fracRelevantReads = reads;
dataIn = struct;

for i=1:6
  tmp = load([bcs_data_path,'MichaelLane/data/sample_',num2str(sampleNum),'_region_',num2str(i),'_unireads.mat']);
  reads{i}.uniqueReads = pack_seqs(tmp.readsuni,64);
  reads{i}.uniqueReads_length = tmp.frequni;
  
  [fracRelevantReads_f{i},sumRelevantReads_f{i}] = currReadsFourth(reads{i}.uniqueReads,reads{i}.uniqueReads_length,values_f{i},0,dataIn);

  [fracRelevantReads_r{i},sumRelevantReads_r{i}] = currReadsFourth(reads{i}.uniqueReads,reads{i}.uniqueReads_length,values_r{i},0,dataIn);
  
end

for i=1:6  
  dat_f(i) = sumRelevantReads_f{i}./sum(reads{i}.uniqueReads_length);
  dat_r(i) = sumRelevantReads_r{i}./sum(reads{i}.uniqueReads_length);
end
%return

dat_f = dat_f';
dat_r = dat_r';
