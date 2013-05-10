function [s,sum_f,sum_r,num_bact_like_W_f,num_bact_like_W_r]=checkForWolbachia(sampleNum,dataBaseFile,w)
%keyboard
bcs_data_path = '~/CS/BAC/';

load(dataBaseFile); 

clear fracRelevantReads sumRelevantReads reads F tmpMat
reads = cell(6,1);
fracRelevantReads = reads;
dataIn = struct;


readLength = 76;
keyboard
for i=1:6
  tmp = load([bcs_data_path,'MichaelLane/data/sample_',num2str(sampleNum),'_region_',num2str(i),'_unireads.mat']);
  
  [junk,i1_f,i2_f] = intersect(tmp.readsuni,seq_readLength76(w,(i-1)*readLength+1:i*readLength),'rows');
  [junk,i1_r,i2_r] = intersect(tmp.readsuni,seq_readLength76(w,(i-1)*readLength+1+readLength*6:i*readLength+readLength*6),'rows');
  
  sum_f(i) = tmp.frequni(i1_f);
  sum_r(i) = tmp.frequni(i1_r);

  s(i) = (sum_f(i)+sum_r(i))/sum(tmp.frequni);
end

%keyboard
% check how many other bacteria match the Wolbachia in each region
for i=1:6 
  sr = sortrows(seq_readLength76(:,(i-1)*readLength+1:i*readLength));
  
  other_f = repmat(seq_readLength76(w,(i-1)*readLength+1:i*readLength),size(sr,1),1);
  
  num_bact_like_W_f(i) = length(find(sum(abs(sr-other_f),2)==0));

  sr = sortrows(seq_readLength76(:,(i-1)*readLength+1+readLength*6:i*readLength+readLength*6));
  other_r = repmat(seq_readLength76(w,(i-1)*readLength+1+readLength*6:i*readLength+readLength*6),size(sr,1),1);
  
  num_bact_like_W_r(i) = length(find(sum(abs(sr-other_r),2)==0));
end

