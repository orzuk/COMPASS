
% check rank in each region

clear
%%%%%%%%%%%55
bcs_data_path = '~/CS/BAC/'
load([bcs_data_path,'full16S/data_primersAmit_2MM_databaseOf2012_readLength_76_noACGT'])

seq_database_file = fullfile(bcs_data_path, '/full16S/data_primersAmit_2MM_databaseOf2012_readLength_76_noACGT_onlySS'); % Database with all sequences
load(seq_database_file); 


sampleNum = 1;

clear fracRelevantReads sumRelevantReads reads F tmpMat
reads = cell(6,1);
fracRelevantReads = reads;
dataIn = struct;
for i=1:6
  tmp = load([bcs_data_path,'MichaelLane/data/sample_',num2str(sampleNum),'_region_',num2str(i),'_unireads.mat']);
  reads{i}.uniqueReads = pack_seqs(tmp.readsuni,64);
  reads{i}.uniqueReads_length = tmp.frequni;
  
  [fracRelevantReads{i},sumRelevantReads{i}] = currReadsFourth(reads{i}.uniqueReads,reads{i}.uniqueReads_length,values{i},0,dataIn);

  F{i} = fracRelevantReads{i}./sum(fracRelevantReads{i});
  
  tmpMat{i} = sparse(size(fracRelevantReads{i},1),6);
  tmpMat{i}(:,i) = -F{i};
end

sumRelevantReads

for i=1:6
  [i,sumRelevantReads{i}./sum(reads{i}.uniqueReads_length)]
  
end



[res,correct,cost,res2] = solve6Regions(M,values,tmpMat,F,SS.Header_uni);


clear read*
sampleNum = 1;
for i=1:6
  tmp = load([bcs_data_path,'MichaelLane/data/sample_',num2str(sampleNum),'_region_',num2str(i),'_unireads.mat']);
  reads8{i}.uniqueReads = pack_seqs(tmp.readsuni,64);
  reads8{i}.uniqueReads_length = tmp.frequni;
  
  [fracRelevantReads8{i},sumRelevantReads8{i}] = currReadsFourth(reads8{i}.uniqueReads,reads8{i}.uniqueReads_length,values{i},0,dataIn);
  
  F8{i} = fracRelevantReads8{i}./sum(fracRelevantReads8{i});
end

sampleNum = 2;
for i=1:6
  tmp = load([bcs_data_path,'MichaelLane/data/sample_',num2str(sampleNum),'_region_',num2str(i),'_unireads.mat']);
  reads16{i}.uniqueReads = pack_seqs(tmp.readsuni,64);
  reads16{i}.uniqueReads_length = tmp.frequni;
  
  [fracRelevantReads16{i},sumRelevantReads16{i}] = currReadsFourth(reads16{i}.uniqueReads,reads16{i}.uniqueReads_length,values{i},0,dataIn);
  F16{i} = fracRelevantReads16{i}./sum(fracRelevantReads16{i});

end

for i=1:6
  [i,sumRelevantReads16{i}./sum(reads16{i}.uniqueReads_length),sumRelevantReads8{i}./sum(reads8{i}.uniqueReads_length)]
  
end




clear s ts
for i=1:6
  [junk,i1,i2] = intersect(reads8{i}.uniqueReads,reads16{i}.uniqueReads,'rows');
  s(i,:) = [sum(reads8{i}.uniqueReads_length(i1)) sum(reads16{i}.uniqueReads_length(i2))];
  ts(i,:) = [sum(reads8{i}.uniqueReads_length) sum(reads16{i}.uniqueReads_length)];
end


for i=1:6
  a = find(F{i}==max(F{i}));
  find(M{i}(a,:))
end

%%%%%%%%%%%%%%%%%%%%
% second try:

clear
bcs_data_path = '~/CS/BAC/';
dataBaseFile = [bcs_data_path,'full16S/data_primersAmit_2MM_databaseOf2012_readLength_76_noACGT'];
seq_database_file = fullfile(bcs_data_path, '/full16S/data_primersAmit_2MM_databaseOf2012_readLength_76_noACGT_onlySS'); % Database with all sequences
load(seq_database_file); 

for sampleNum=1:4
  sampleNum
  solveMichaelFirstTry(sampleNum,dataBaseFile);
end

for sampleNum=[1:4]
  sampleNum
  x = load([bcs_data_path,'MichaelLane/results/firstTryMichael_',num2str(sampleNum)]);
  fileName = [bcs_data_path,'MichaelLane/results/res_firstTryMichael_',num2str(sampleNum),'.xls'];
  printResSeq(x.b,x.results,SS.Header_uni,fileName);
  fclose('all')
  
  %pause
end

% Wolbachia
for i=1:length(SS.Header_uni)
  for j=1:length(SS.Header_uni{i})
    if iscell(SS.Header_uni{i}{j})
      for k=1:length(SS.Header_uni{i}{j})
        f = findstr(SS.Header_uni{i}{j}{k},'108447 NC_002978.6');
        if ~isempty(f)
          SS.Header_uni{i}{j}{k}
          [i]
          pause
        end
      end
    else
      f = findstr(SS.Header_uni{i}{j},'108447 NC_002978.6');
      if ~isempty(f)
        SS.Header_uni{i}{j}
        [i]
        pause
      end
    end
    
    
  end
end




% check for Wolbachia
w  = 542133;
clear s sum_f sum_r num_bact_like_W_f num_bact_like_W_r
for sampleNum=1:4
  [s{sampleNum},sum_f{sampleNum},sum_r{sampleNum},num_bact_like_W_f{sampleNum},num_bact_like_W_r{sampleNum}] = checkForSpecificBacteria(sampleNum,dataBaseFile,w,M,values,seq_readLength76,readLength);
end


%%%%%%%%%%%%


%%%%
% Lactobacillus

nm_454 = {x...
'282170 AB368903.1',...
'73879 AJ272031.1' ,...
'560153 FJ455519.1',... 
'161334 AM279762.2' ,...
'542366 GQ922601.2'};

ind_454 = cell(1,length(nm_454));
for i=1:length(nm_454)
  ind_454{i} = findSpecificBacteria2(SS.Header_uni,nm_454{i});
end

readLength = 76;
load(dataBaseFile)

clear s_454
sum_f_454 sum_r_454 num_bact_like_W_f_454 num_bact_like_W_r_454
for sampleNum=1:4
  for i=1:length(nm_454)
    [s_454{sampleNum}(i,:),sum_f_454{sampleNum}(i,:),sum_r_454{sampleNum}(i,:),num_bact_like_W_f_454{sampleNum}(i,:),num_bact_like_W_r_454{sampleNum}(i,:)] = checkForSpecificBacteria(sampleNum,dataBaseFile,ind_454{i},M,values,seq_readLength76,readLength);
  end
end


s_454

sum_f_454

sum_r_454

num_bact_like_W_f_454

num_bact_like_W_r_454


nm_Illumina = {...
'546229 GU195646.1' ,...
'284270 AB368910.1',... 
'587695 GU138567.1',... 
'106977 AJ640078.1',... 
'161334 AM279762.2',... 
'583941 GU138597.1',... 
'534907 GQ461598.1',... 
'290616 EU675926.1',... 
'15163 M58827.1' ,...
'548471 GQ167190.1' ,...
'257487 EU081013.1' ,...
'551236 GQ202837.1'};


for sampleNum=1:4
  sampleNum
  solveMichaelFirstTry_partOfRegions(sampleNum,dataBaseFile);
end





%%%%%%%%%%%%%%%%5
%%%%
% third try
clear
bcs_data_path = '~/CS/BAC/';
dataBaseFile = [bcs_data_path,'full16S/data_primersAmit_0MM_databaseOf2012_readLength_76_noACGT'];
load(dataBaseFile)
readLength = 76;
seq_database_file = fullfile(bcs_data_path, '/full16S/data_primersAmit_0MM_databaseOf2012_readLength_76_noACGT_onlySS'); % Database with all sequences
load(seq_database_file); 



for sampleNum=5:20
  sampleNum
  solveMichaelThirdTry(sampleNum,dataBaseFile);
end

for sampleNum=[5:20]
  sampleNum
  x = load([bcs_data_path,'MichaelLane/results/thirdTryMichael_',num2str(sampleNum)]);
  fileName = [bcs_data_path,'MichaelLane/results/res_thirdTryMichael_',num2str(sampleNum),'.xls'];
  printResSeq(x.b,x.results,SS.Header_uni,fileName);
  fclose('all')
  
  %pause
end

% Wolbachia
for i=1:length(SS.Header_uni)
  for j=1:length(SS.Header_uni{i})
    if iscell(SS.Header_uni{i}{j})
      for k=1:length(SS.Header_uni{i}{j})
        f = findstr(SS.Header_uni{i}{j}{k},'108447 NC_002978.6');
        if ~isempty(f)
          SS.Header_uni{i}{j}{k}
          [i]
          pause
        end
      end
    else
      f = findstr(SS.Header_uni{i}{j},'108447 NC_002978.6');
      if ~isempty(f)
        SS.Header_uni{i}{j}
        [i]
        pause
      end
    end
    
    
  end
end

% check for Wolbachia
w  = 223689;
clear s sum_f sum_r num_bact_like_W_f num_bact_like_W_r
for sampleNum=1:4
  [s{sampleNum},sum_f{sampleNum},sum_r{sampleNum},num_bact_like_W_f{sampleNum},num_bact_like_W_r{sampleNum}] = checkForSpecificBacteria(sampleNum,dataBaseFile,w,M,values,seq_readLength76,readLength);
end

% look at results of all samples - see reoccuring bacteria
% find concensus of reads that appear in all cases - are they mapped?
% how differnet are the reads that were not mapped

[cscore,algn]=swalign(seq_tmpInd,seq_rec,'Alphabet','NT');

%%%%%%%%%%%%5
clear
bcs_data_path = '~/CS/BAC/';
load([bcs_data_path,'full16S/data_primersAmit_0MM_databaseOf2012_readLength_76_noACGT_M_values_forward_reverse_separated'])

dataBaseFile = [bcs_data_path,'full16S/data_primersAmit_0MM_databaseOf2012_readLength_76_noACGT_M_values_forward_reverse_separated'];

for sampleNum=1:20
  [dat_f(sampleNum,:),dat_r(sampleNum,:)] = compareMichael_forward_and_reverse(sampleNum,dataBaseFile);
end

for sampleNum=1:20
  [dat(sampleNum,:)] = compareMichael_takenTogether(sampleNum,[bcs_data_path,'full16S/data_primersAmit_0MM_databaseOf2012_readLength_76_noACGT']);
end



figure(1)
for i=1:6
  subplot(3,2,i);
  plot(dat_f(:,i),dat_r(:,i),'.');
  xlabel('forward');ylabel('reverse')
  title(i)
  set(gca,'xlim',[0 1],'ylim',[0 1])
  line([0 0.5],[0 0.5])
end

print -dpdf ~/CS/BAC/MichaelLane/results/alignmentPercentage

figure
plot(dat_f+dat_r,dat,'.')

figure(2)
for i=1:20
  subplot(4,5,i);
  plot(dat_f(i,:),dat_r(i,:),'.');
  title(i)
  set(gca,'xlim',[0 1],'ylim',[0 1])
  line([0 0.5],[0 0.5])
end


[x,y] = find(dat<0.3 & dat_f+dat_r>0.9);



sampleNum = 9

i = 6;
tmp = load([bcs_data_path,'MichaelLane/data/sample_',num2str(sampleNum),'_region_',num2str(i),'_unireads.mat']);
reads{i}.uniqueReads = pack_seqs(tmp.readsuni,64);
reads{i}.uniqueReads_length = tmp.frequni;

c = load([bcs_data_path,'full16S/data_primersAmit_0MM_databaseOf2012_readLength_76_noACGT']);
fr = load([bcs_data_path,'full16S/data_primersAmit_0MM_databaseOf2012_readLength_76_noACGT_M_values_forward_reverse_separated']);

dataIn = struct;
[fracRelevantReads_c{i},sumRelevantReads_c{i}] = currReadsFourth(reads{i}.uniqueReads,reads{i}.uniqueReads_length,c.values{i},0,dataIn);

[fracRelevantReads_f{i},sumRelevantReads_f{i}] = currReadsFourth(reads{i}.uniqueReads,reads{i}.uniqueReads_length,fr.values_f{i},0,dataIn);

[fracRelevantReads_r{i},sumRelevantReads_r{i}] = currReadsFourth(reads{i}.uniqueReads,reads{i}.uniqueReads_length,fr.values_r{i},0,dataIn);


%%%%%%%%5
clear
bcs_data_path = '~/CS/BAC/';
dataBaseFile = [bcs_data_path,'full16S/data_primersAmit_0MM_databaseOf2012_readLength_76_noACGT'];
%load(dataBaseFile)
%readLength = 76;
seq_database_file = fullfile(bcs_data_path, '/full16S/data_primersAmit_0MM_databaseOf2012_readLength_76_noACGT_onlySS'); % Database with all sequences
load(seq_database_file); 



for sampleNum=1:20
  sampleNum
  %solveMichaelThirdTry(sampleNum,dataBaseFile);
  solveMichaelFourthTry(sampleNum,dataBaseFile);
end

for sampleNum=1:20
  sampleNum
  x = load([bcs_data_path,'MichaelLane/results/thirdTryMichael_',num2str(sampleNum)]);
  fileName = [bcs_data_path,'MichaelLane/results/res_thirdTryMichael_',num2str(sampleNum),'.xls'];
  printResSeq(x.b,x.results,SS.Header_uni,fileName);
  fclose('all')
  
  %pause
end

%%%%%%%%%%%%%%%%
% add 1 random mixtures
clear 
bcs_data_path = '~/CS/BAC/'
dataBaseFile = [bcs_data_path,'full16S/data_primersAmit_0MM_databaseOf2012_readLength_76_noACGT'];
load(dataBaseFile)

readLength = 76;

a = randperm(size(seq_readLength76,1));
sel = a(1:10);
sampleNum = 21;
for i=1:6
    reg = [];
    for j=sel
      if seq_readLength76(j,(i-1)*readLength+1)~=char(1)
        reg = [reg;seq_readLength76(j,(i-1)*readLength+1:i*readLength);seq_readLength76(j,(i-1)*readLength+1+readLength*6:i*readLength+readLength*6)];
      end
      
    end
    
    red = pack_seqs(reg,64);

    [uniqueReads,uniqueReads_inds] = extract_sub_kmers(red, readLength*ones(size(red,1),1),readLength, 1,0);
    clear red
    
    [junk_vals junk_inds uniqueReads_length] = get_duplicates(uniqueReads_inds(:,1));
    clear junk_vals junk_inds uniqueReads_inds
    
    readsuni = int2nt(unpack_seqs(uniqueReads,readLength*ones(1,size(uniqueReads,1)),64));
    frequni = uniqueReads_length;
    save([bcs_data_path,'MichaelLane/data/sample_',num2str(sampleNum),'_region_',num2str(i),'_unireads.mat'],'readsuni','frequni');
end

solveMichaelFourthTry(21,dataBaseFile);


% find Wolbachia and bumble bee 

seq_database_file = fullfile(bcs_data_path, '/full16S/data_primersAmit_0MM_databaseOf2012_readLength_76_noACGT_onlySS'); % Database with all sequences
load(seq_database_file); 


findSpecificBacteria2(SS.Header_uni,'827815 HM109722.1')
bb_0MM = 18951

s_bb = seq_readLength76(bb_0MM,:);

findSpecificBacteria2(SS.Header_uni,'273831 EU096232.1')
wol_0MM = 223689

s_wol = seq_readLength76(wol_0MM,:);


save ~/CS/BAC/MichaelLane/results/checkSeq s_bb s_wol


%%%%%%%%%%%%%
% misclassification of reginos

clear
bcs_data_path = '~/CS/BAC/';
dataBaseFile = [bcs_data_path,'full16S/data_primersAmit_0MM_databaseOf2012_readLength_76_noACGT'];
load(dataBaseFile); 

clear fracRelevantReads sumRelevantReads reads F tmpMat
reads = cell(6,1);
fracRelevantReads = reads;
dataIn = struct;

sampleNum = 5;

i = 2;
tmp = load([bcs_data_path,'MichaelLane/data/sample_',num2str(sampleNum),'_region_',num2str(i),'_unireads.mat']);
reads{i}.uniqueReads = pack_seqs(tmp.readsuni,64);
reads{i}.uniqueReads_length = tmp.frequni;

for i=1:6
  [f,sumRelevantReads] ...
      = currReadsFourth(reads{2}.uniqueReads,reads{2}.uniqueReads_length,values{i},0,dataIn);
  sumRelevantReads
end

%%%%%%%%%%%
clear
bcs_data_path = '~/CS/BAC/';
dataBaseFile = [bcs_data_path,'full16S/data_primersAmit_0MM_databaseOf2012_readLength_76_noACGT'];
seq_database_file = fullfile(bcs_data_path, '/full16S/data_primersAmit_0MM_databaseOf2012_readLength_76_noACGT_onlySS'); % Database with all sequences
load(seq_database_file); 

for sampleNum=[5:15]
  sampleNum
  dat(sampleNum,:)=solveMichaelFourthTry(sampleNum,dataBaseFile);
end

for sampleNum=13]
  sampleNum
  x = load([bcs_data_path,'MichaelLane/results/fourthTryMichael_',num2str(sampleNum)]);
  fileName = [bcs_data_path,'MichaelLane/results/res_fourthTryMichael_',num2str(sampleNum),'.xls'];
  printResSeq(x.b,x.results,SS.Header_uni,fileName);
  fclose('all')
  
  %pause
end


%%%%%%% align issues:
clear
bcs_data_path = '~/CS/BAC/';
dataBaseFile = [bcs_data_path,'full16S/data_primersAmit_0MM_databaseOf2012_readLength_76_noACGT'];

for sampleNum=[1:20]
  sampleNum
  dat0(sampleNum,:)=solveMichaelFourthTry(sampleNum,dataBaseFile);
end

figure(1);imagesc(dat0);title('0MM');colorbar
print -dpdf ~/CS/BAC/MichaelLane/results/alignment_0MM

bcs_data_path = '~/CS/BAC/';
dataBaseFile = [bcs_data_path,'full16S/data_primersAmit_2MM_databaseOf2012_readLength_76_noACGT'];
for sampleNum=[1:20]
  sampleNum
  dat2(sampleNum,:)=solveMichaelFourthTry(sampleNum,dataBaseFile);
end
figure(2);imagesc(dat2);title('2MM');colorbar
print -dpdf ~/CS/BAC/MichaelLane/results/alignment_2MM



findSpecificBacteria2(SS.Header_uni,'NC_007482')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
bcs_data_path = '~/CS/BAC/';
dataBaseFile = [bcs_data_path,'full16S/data_primersAmit_0MM_databaseOf2012_readLength_76_noACGT'];

for sampleNum=[1:20]
  sampleNum
  solveMichaelFourthTry(sampleNum,dataBaseFile);
end
