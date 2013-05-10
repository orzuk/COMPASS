%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read length 75
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
bcs_data_path = '~/CS/BAC/';
seq_database_file = fullfile(bcs_data_path, '/full16S/bac16s_full_without_ambiguous'); % Database with all sequences
load([bcs_data_path,'full16S/primerDataAmitPrimers_full_without_ambiguous'])

readLength = 75; % set read length  - this corresponds to 100 !!!!
AmitBarcodeOffset = 18; % start read from this.


primer_fields = fields(P_no_nonACGT); 
seq_ind = zeros(size(P_no_nonACGT.pos_pr1,1),length(primer_fields));
pos_f = seq_ind;
pos_r = seq_ind;
for i=1:length(primer_fields)
  
  w = ['amp = find(P_no_nonACGT.',primer_fields{i},'(:,1)~=0 & P_no_nonACGT.',primer_fields{i},'(:,2)~=0);'];
  eval(w);
  seq_ind(amp,i) = amp;

  w = ['f = P_no_nonACGT.' primer_fields{i} '(amp,1)+AmitBarcodeOffset;'];
  eval(w);
  pos_f(amp,i) = f;
  
  w = ['r = P_no_nonACGT.' primer_fields{i} '(amp,2)-((readLength+AmitBarcodeOffset-1));'];
  eval(w);
  pos_r(amp,i) = r;
  
end
 

%%%%%%%%%%
if 1==1 % prepare the sequences
  S = load(seq_database_file); 
  S.Sequence_uni = vec2column(S.Sequence_uni); 
  
  % create concatenated sequence
  seq_f = char(ones(size(pos_f,1),6*readLength));
  seq_r = seq_f;
  
  for i=1:length(pos_f)
    if mod(i,1000)==1
      i
    end
    a = find(pos_f(i,:));
    for j=a
      seq_f(i,(j-1)*readLength+1:j*readLength) = S.Sequence_uni{i}(pos_f(i,j):pos_f(i,j)+readLength-1);
      seq_r(i,(j-1)*readLength+1:j*readLength) = S.Sequence_uni{i}(pos_r(i,j):pos_r(i,j)+readLength-1);
    end
  end
  
  [tmp_seq_readLength75,origInd] = sortrows([seq_f,seq_r]);


  [u_f,i1_f] = unique(tmp_seq_readLength75,'rows','first');
  [u_l,i1_l] = unique(tmp_seq_readLength75,'rows','last');
  
  seq_readLength75 = u_f;
  
  % find the bacteria that are part of this 
  sameSeq = cell(length(i1_f),1);
  group = zeros(length(S.Sequence_uni),1);
  for i=1:length(i1_f)
    sameSeq{i} = origInd(i1_f(i):i1_l(i));
    group(sameSeq{i}) = i;
  end
  
  
  
end % create the reads 



% matrix 

% prepare matrix for each region
load  ~/CS/BAC/full16S/Sequence_packed_full_without_ambiguous

% check whether each region is unique

clear M value
[group_val,inds]= unique(group);  % take a representative
for reg=1:6
  % output the bacteria in each region
  [PrimerReadBySpeciesMat_tmp Primer_kmers_packed_tmp bactInRegion_tmp]=oneRegionData(reg,inds,seq_ind,pos_f,pos_r,readLength,SS);
  
  [junk,i1] = intersect(inds,bactInRegion_tmp);
  
  M{reg} = sparse(size(PrimerReadBySpeciesMat_tmp,1),length(inds));
  M{reg}(:,i1) = PrimerReadBySpeciesMat_tmp(:,bactInRegion_tmp);
  values{reg} = Primer_kmers_packed_tmp;
end

save([bcs_data_path,'full16S/data_primersAmit_0MM_readLength_75'],'sameSeq','seq_readLength75','group','M','values','inds') 

%check - overlap bwtwqeen region
for i=1:5
  for j=i+1:6
    junk=intersect(values{i},values{j},'rows');
    if ~isempty(junk)
      disp('overlap!!!')
    end
  end
end


if 1==2
  bcs_data_path = '~/CS/BAC/';
  load([bcs_data_path,'full16S/data_primersAmit_0MM_readLength_75'])

  
  [x,y] = find(seq_readLength75==char(1));
  rest = 1:size(seq_readLength75,1);
  rest(x) = [];
  
  less = 1:size(seq_readLength75,1);
  less(rest) = [];
  
  clear reads
  numBact = 5;
  for i=1:6
    reg = [];
    for j=rest(1:numBact)
      if seq_readLength75(j,(i-1)*75+1)~=char(1)
        reg = [reg;seq_readLength75(j,(i-1)*75+1:i*75);seq_readLength75(j,(i-1)*75+451:i*75+450)];
      end
      
    end
    for j=less(10)
      if seq_readLength75(j,(i-1)*75+1)~=char(1)
        reg = [reg;repmat([seq_readLength75(j,(i-1)*75+1:i*75);seq_readLength75(j,(i-1)*75+451:i*75+450)],10,1)];
      end
      
    end
    red = pack_seqs(reg,64);

    [uniqueReads,uniqueReads_inds] = extract_sub_kmers(red, 75*ones(size(red,1),1),75, 1,0);
    clear red
    
    [junk_vals junk_inds uniqueReads_length] = get_duplicates(uniqueReads_inds(:,1));
    clear junk_vals junk_inds uniqueReads_inds
    
    
    reads{i}.uniqueReads = uniqueReads;
    reads{i}.uniqueReads_length = uniqueReads_length;
  end
  
  
  
  clear reads
  numBact = 10;
  for i=1:6
    
    red = pack_seqs([seq_readLength75(rest(1:numBact),(i-1)*75+1:i*75);seq_readLength75(rest(1:numBact),(i-1)*75+451:i*75+450)],64);
    
    red = [red;pack_seqs(repmat([seq_readLength75(rest(numBact),(i-1)*75+1:i*75);seq_readLength75(rest(numBact),(i-1)*75+451:i*75+450)],5,1),64)];
    [uniqueReads,uniqueReads_inds] = extract_sub_kmers(red, 75*ones(size(red,1),1),75, 1,0);
    clear red
    
    [junk_vals junk_inds uniqueReads_length] = get_duplicates(uniqueReads_inds(:,1));
    clear junk_vals junk_inds uniqueReads_inds
    
    
    reads{i}.uniqueReads = uniqueReads;
    reads{i}.uniqueReads_length = uniqueReads_length;
  end
  
  clear fracRelevantReads tmpMat F
  for i=1:6
    [fracRelevantReads{i},sumRelevantReads{i}] = currReadsFourth(reads{i}.uniqueReads,reads{i}.uniqueReads_length,values{i},0,struct);
    
    F{i} = fracRelevantReads{i}./sum(fracRelevantReads{i});
    
    tmpMat{i} = sparse(size(fracRelevantReads{i},1),6);
    tmpMat{i}(:,i) = -F{i};
  end

  
  Ap=[M{1},tmpMat{1};M{2},tmpMat{2};M{3},tmpMat{3};M{4},tmpMat{4};M{5},tmpMat{5};M{6},tmpMat{6}];
  [x] = solve_L2_2(Ap,zeros(size(Ap,1),1));

  res = x(1:end-6);
  res = res./sum(res);

  find(res)
  
end


clear res
for i=1:6
  
  Ap=[M{i},tmpMat{i}];
  [x] = solve_L2_2(Ap,zeros(size(Ap,1),1));

  tmp_res = x(1:end-1);
  res{i} = tmp_res./sum(tmp_res);

  %find(res)

end




clear reads
numBact = 2;

for i=1:6
        red = pack_seqs([seq_readLength75(rest(1:numBact),(i-1)*75+1:i*75);seq_readLength75(rest(1:numBact),(i-1)*75+451:i*75+450)],64);
    
    [uniqueReads,uniqueReads_inds] = extract_sub_kmers(red, 75*ones(size(red,1),1),75, 1,0);
    clear red
    
    [junk_vals junk_inds uniqueReads_length] = get_duplicates(uniqueReads_inds(:,1));
    clear junk_vals junk_inds uniqueReads_inds
    
    
    reads{i}.uniqueReads = uniqueReads;
    reads{i}.uniqueReads_length = uniqueReads_length+0.01*randn(size(uniqueReads_length));
end 

for i=1:6
  
    [fracRelevantReads{i},sumRelevantReads{i}] = currReadsFourth(reads{i}.uniqueReads,reads{i}.uniqueReads_length,values{i},0,struct);
    
    F{i} = fracRelevantReads{i}./sum(fracRelevantReads{i});
    
    tmpMat{i} = sparse(size(fracRelevantReads{i},1),6);
    tmpMat{i}(:,i) = -F{i};
end
Ap=[M{1},tmpMat{1};M{2},tmpMat{2};M{3},tmpMat{3};M{4},tmpMat{4};M{5},tmpMat{5};M{6},tmpMat{6}];
[x] = solve_L2_2(Ap,zeros(size(Ap,1),1));

res = x(1:end-6);
res = res./sum(res);


%%%%%%%%%%%%%%%%%
a = randperm(size(seq_readLength75,1));
sel = a(1:10);
for i=1:6
    reg = [];
    for j=sel
      if seq_readLength75(j,(i-1)*75+1)~=char(1)
        reg = [reg;seq_readLength75(j,(i-1)*75+1:i*75);seq_readLength75(j,(i-1)*75+451:i*75+450)];
      end
      
    end
    
    red = pack_seqs(reg,64);

    [uniqueReads,uniqueReads_inds] = extract_sub_kmers(red, 75*ones(size(red,1),1),75, 1,0);
    clear red
    
    [junk_vals junk_inds uniqueReads_length] = get_duplicates(uniqueReads_inds(:,1));
    clear junk_vals junk_inds uniqueReads_inds
    
    
    reads{i}.uniqueReads = uniqueReads;
    reads{i}.uniqueReads_length = uniqueReads_length;
end

for i=1:6
  
    [fracRelevantReads{i},sumRelevantReads{i}] = currReadsFourth(reads{i}.uniqueReads,reads{i}.uniqueReads_length,values{i},0,struct);
    
    F{i} = fracRelevantReads{i}./sum(fracRelevantReads{i});
    
    tmpMat{i} = sparse(size(fracRelevantReads{i},1),6);
    tmpMat{i}(:,i) = -F{i};
end

Ap=[M{1},tmpMat{1};M{2},tmpMat{2};M{3},tmpMat{3};M{4},tmpMat{4};M{5},tmpMat{5};M{6},tmpMat{6}];
[x] = solve_L2_2(Ap,zeros(size(Ap,1),1));

res = x(1:end-6);
res = res./sum(res);


%%%%%%%%%%

%sel = [136065 185801 85093 200836 98816 28845 15557 31485 19262 65477];
%sel = [124258 63421       98470      196140       69392      181451       81548      172128      169865     76630];
sel = [36065  85093       172128]
for i=1:6
    reg = [];
    for j=sel
      if seq_readLength75(j,(i-1)*75+1)~=char(1)
        reg = [reg;seq_readLength75(j,(i-1)*75+1:i*75);seq_readLength75(j,(i-1)*75+451:i*75+450)];
      else
        %no = [no;i,j];
      end
      
    end
    
    red = pack_seqs(reg,64);

    [uniqueReads,uniqueReads_inds] = extract_sub_kmers(red, 75*ones(size(red,1),1),75, 1,0);
    clear red
    
    [junk_vals junk_inds uniqueReads_length] = get_duplicates(uniqueReads_inds(:,1));
    clear junk_vals junk_inds uniqueReads_inds
    
    
    reads{i}.uniqueReads = uniqueReads;
    reads{i}.uniqueReads_length = uniqueReads_length;
end


for i=1:6
    [fracRelevantReads{i},sumRelevantReads{i}] = currReadsFourth(reads{i}.uniqueReads,reads{i}.uniqueReads_length,values{i},0,struct);
    
    F{i} = fracRelevantReads{i}./sum(fracRelevantReads{i});
    
    tmpMat{i} = sparse(size(fracRelevantReads{i},1),6);
    tmpMat{i}(:,i) = -F{i};
end

%Ap=[M{1},tmpMat{1};M{2},tmpMat{2};M{3},tmpMat{3};M{4},tmpMat{4};M{5},tmpMat{5};M{6},tmpMat{6}];
%[x] = solve_L2_2(Ap,zeros(size(Ap,1),1));

%Ap=[M{1},tmpMat{1};M{2},tmpMat{2};M{3},tmpMat{3};M{4},tmpMat{4};M{5},tmpMat{5};M{6},tmpMat{6}];
%Ap=[Ap;10^-4*ones(1,size(M{1},2),1),zeros(1,6)];

%[x] = solve_L2_2(Ap,[zeros(size(Ap,1)-1,1);10^-4]);


%Ap=[M{1},tmpMat{1};M{2},tmpMat{2};M{3},tmpMat{3};M{4},tmpMat{4};M{5},tmpMat{5};M{6},tmpMat{6}];
%Ap=[ones(1,size(M{1},2),1),zeros(1,6);Ap];
%[x] = solve_L2_2(Ap,[1;zeros(size(Ap,1)-1,1)]);


%Ap=[M{1},tmpMat{1};M{2},tmpMat{2};M{3},tmpMat{3};M{4},tmpMat{4};M{5},tmpMat{5};M{6},tmpMat{6}];

Ap(:,1) = []; % do we need this?
[x] = solve_L2_3(Ap,zeros(size(Ap,1),1));

res = x(1:end-6);

%%%%%%%%%%%%%%%%
clear

bcs_data_path = '~/CS/BAC/';
load([bcs_data_path,'full16S/data_primersAmit_0MM_readLength_75'])

load ~/CS/BAC/test6Region/sim_6regions_bac_dist_flag_1_Nbac_in_mixture_200_npower_100_readlen_75_numIter_1_Nreads_100000_Noise_0_Correction_0

% test infinite reads - with pruning

%sel = [36065  85093       172128]
sel = [36065       172128]
[junk i1 i2] = intersect(bact,sel);

%specificBAC6Regions(bact(i1),freq(i1),M,seq_readLength75,values)
      

sel = unique([36065       172128      ]);
[res,correct,cost,res2]=specificBAC6Regions(sel,[0.5 0.5],M,seq_readLength75,values);

[res,correct,cost,res2]=specificBAC6Regions_Sec(sel,[0.8 0.2],M,seq_readLength75,values);

[res,correct,cost,res2]=specificBAC6Regions_Third(sel,[0.8 0.2],M,seq_readLength75,values);

specificBAC6Regions(bact,freq,M,seq_readLength75,values)

sm = [36065       172128 bact(2:5)];bact(2)
%specificBAC6Regions_Sec(sm,1/length(sm)*ones(1,length(sm)),M,seq_readLength75,values)


sm8 = ceil(rand(1,100)*length(sameSeq));
reads8 = specificBAC6Regions_junk(sm8,1/length(sm8)*ones(1,length(sm8)),M,seq_readLength75,values)
sm16 = ceil(rand(1,100)*length(sameSeq));
reads16 = specificBAC6Regions_junk(sm16,1/length(sm16)*ones(1,length(sm16)),M,seq_readLength75,values)

clear s ts
for i=1:6
  [junk,i1,i2] = intersect(reads8{i}.uniqueReads,reads16{i}.uniqueReads,'rows');
  s(i,:) = [sum(reads8{i}.uniqueReads_length(i1)) sum(reads16{i}.uniqueReads_length(i2))];
  ts(i,:) = [sum(reads8{i}.uniqueReads_length) sum(reads16{i}.uniqueReads_length)];
end

s./ts

keyboard

