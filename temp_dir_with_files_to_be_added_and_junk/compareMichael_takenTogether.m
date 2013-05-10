function dat=compareMichael_takenTogether(sampleNum,dataBaseFile)
%keyboard
bcs_data_path = '~/CS/BAC/';

load(dataBaseFile); 
%load('tmp')

clear fracRelevantReads sumRelevantReads reads F tmpMat
reads = cell(6,1);
fracRelevantReads = reads;
dataIn = struct;

for i=1:6
  tmp = load([bcs_data_path,'MichaelLane/data/sample_',num2str(sampleNum),'_region_',num2str(i),'_unireads.mat']);
  reads{i}.uniqueReads = pack_seqs(tmp.readsuni,64);
  reads{i}.uniqueReads_length = tmp.frequni;
  
  [fracRelevantReads{i},sumRelevantReads{i}] = currReadsFourth(reads{i}.uniqueReads,reads{i}.uniqueReads_length,values{i},0,dataIn);

end

for i=1:6  
  dat(i) = sumRelevantReads{i}./sum(reads{i}.uniqueReads_length);
end

dat = dat';