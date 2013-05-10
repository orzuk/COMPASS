
% check rank in each region

clear
%%%%%%%%%%%55
bcs_data_path = '~/CS/BAC/'
load([bcs_data_path,'full16S/data_primersAmit_2MM_databaseOf2012_readLength_76_noACGT_forPE'])

seq_database_file = fullfile(bcs_data_path, '/full16S/data_primersAmit_2MM_databaseOf2012_readLength_76_noACGT_onlySS_forPE'); % Database with all sequences
load(seq_database_file); 


sampleNum = 10;

clear fracRelevantReads sumRelevantReads reads F tmpMat
reads = cell(6,1);
fracRelevantReads = reads;
dataIn = struct;
for i=1:6
  tmp = load([bcs_data_path,'Niv/PE76/data/sample_',num2str(sampleNum),'_region_',num2str(i),'_unireads.mat']);
  reads{i}.uniqueReads = pack_seqs(tmp.readsuni,64);
  reads{i}.uniqueReads_length = tmp.frequni;
  
  [fracRelevantReads{i},sumRelevantReads{i}] = currReadsFourth(reads{i}.uniqueReads,reads{i}.uniqueReads_length,values{i},0,dataIn);

  F{i} = fracRelevantReads{i}./sum(fracRelevantReads{i});
  
  tmpMat{i} = sparse(size(fracRelevantReads{i},1),6);
  tmpMat{i}(:,i) = -F{i};
end

for i=1:6  
  [i sumRelevantReads{i}./sum(reads{i}.uniqueReads_length)]
end


%[res,correct,cost,res2] = solve6Regions(M,values,tmpMat,F,SS.Header_uni);

[res,correct,cost,res2] = solveRegion2(M,values,tmpMat,F,SS.Header_uni);




clear read*
sampleNum = 1;
for i=1:6
  tmp = load([bcs_data_path,'Niv/PE76/data/sample_',num2str(sampleNum),'_region_',num2str(i),'_unireads.mat']);
  reads8{i}.uniqueReads = pack_seqs(tmp.readsuni,64);
  reads8{i}.uniqueReads_length = tmp.frequni;
end

sampleNum = 20;
for i=1:6
  tmp = load([bcs_data_path,'Niv/PE76/data/sample_',num2str(sampleNum),'_region_',num2str(i),'_unireads.mat']);
  reads16{i}.uniqueReads = pack_seqs(tmp.readsuni,64);
  reads16{i}.uniqueReads_length = tmp.frequni;
end


clear s ts
for i=1:6
  [junk,i1,i2] = intersect(reads8{i}.uniqueReads,reads16{i}.uniqueReads,'rows');
  s(i,:) = [sum(reads8{i}.uniqueReads_length(i1)) sum(reads16{i}.uniqueReads_length(i2))];
  ts(i,:) = [sum(reads8{i}.uniqueReads_length) sum(reads16{i}.uniqueReads_length)];
end

%%%%%%%%%
% solve and keep the data
clear
bcs_data_path = '~/CS/BAC/';
dataBaseFile = [bcs_data_path,'full16S/data_primersAmit_2MM_databaseOf2012_readLength_76_noACGT_forPE'];

%for sampleNum=[1:16 18:20]
for sampleNum=18:20
  sampleNum
  solveNivFirstTry(sampleNum,dataBaseFile);
end

for sampleNum=[1:16 18:20]
  sampleNum
  x = load([bcs_data_path,'Niv/PE76/results/firstTry_',num2str(sampleNum)]);
  [length(x.b) length(x.b_unif) length(intersect(x.b_unif,x.b)) length(setdiff(x.b_unif,x.b)) length(setdiff(x.b,x.b_unif))]
  x.results_unif
  pause
end
% %%%%%%
