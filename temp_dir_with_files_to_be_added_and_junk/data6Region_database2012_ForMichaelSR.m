%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read length 76
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('never pack with nonACGT !!!!!!')

% create packed sequences

clear
bcs_data_path = '~/CS/BAC/';
seq_database_file = fullfile(bcs_data_path, '/database2012/bac16s_database_uni_download_nov11_2012'); % Database with all sequences
P = load([bcs_data_path,'6region/primers_position_few_regions_2MM_database780000'])

readLength = 76; % set read length  - this corresponds to 100 !!!!
AmitBarcodeOffset = 18; % start read from this.


primer_fields = fields(P); 
seq_ind = zeros(size(P.pos_pr1,1),length(primer_fields));
pos_f = seq_ind;
pos_r = seq_ind;
for i=1:length(primer_fields)
  
  w = ['amp = find(P.',primer_fields{i},'(:,1)~=0 & P.',primer_fields{i},'(:,2)~=0);'];
  eval(w);
  seq_ind(amp,i) = amp;

  w = ['f = P.' primer_fields{i} '(amp,1)+AmitBarcodeOffset;'];
  eval(w);
  pos_f(amp,i) = f;
  
  w = ['r = P.' primer_fields{i} '(amp,2)-((readLength+AmitBarcodeOffset-1));'];
  eval(w);
  pos_r(amp,i) = r;
end
 

%%%%%%%%%%
% prepare the sequences
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

% this includes the nonACGT
[tmp_seq_readLength76_after_sort,origInd] = sortrows([seq_f,seq_r]);

[u_f,i1_f] = unique(tmp_seq_readLength76_after_sort,'rows','first');
[u_l,i1_l] = unique(tmp_seq_readLength76_after_sort,'rows','last');
 
seq_readLength76 = u_f; % still contains nonACGT

% same seq is in the original 780000 database 
sameSeq = cell(size(seq_readLength76,1),1);
for i=1:length(i1_f)
  sameSeq{i} = origInd(i1_f(i):i1_l(i));
end

% delete cases of nonACGT in 76*12
[bact_with_nonACGT_over_76,y] = find(seq_readLength76~='A' & seq_readLength76~='C' & seq_readLength76~='G' & seq_readLength76~='T' & seq_readLength76~=char(1)); 

seq_readLength76(bact_with_nonACGT_over_76,:) = [];
sameSeq(bact_with_nonACGT_over_76) = [];

% the first line is empty
if isempty(find(abs(seq_readLength76(1,:))~=1))
  seq_readLength76(1,:) = [];
  sameSeq(1) = [];
end



% prepare the database for 76 - without the ambiguous over the 76*12

% delete all bact with nonACGT in the readLength*12 region
% create new sequences which have 'A' in the not amplified regions
% pack them
% use the regions to create them

lengthPacked = length(pack_seqs(char('A'*ones(1,readLength*2*6)),64));

numBact = length(sameSeq);

basicOffSet = readLength*6;
clear SS
SS.len_uni = readLength*6*2*ones(1,numBact);
SS.amplifiedRegion = zeros(numBact,6);
SS.Sequence_packed_matrix = zeros(numBact,lengthPacked,'uint64');
SS.Sequence_uni_matrix = char('Z'*zeros(numBact,readLength*6*2));
SS.README = 'all forward and then all reverse';
SS.Header_uni = cell(1,numBact);
curr_time = clock;

for i=1:numBact
  
  if mod(i,1000)==0
    i
    etime(clock,curr_time)
    curr_time = clock;
  end
  
  curr_seq = seq_readLength76(i,:);
  a = find(abs(curr_seq)==1);
  curr_seq(a) = 'A';
  SS.Sequence_uni_matrix(i,:) = curr_seq;
  SS.Sequence_packed_matrix(i,:) = pack_seqs(curr_seq,64);
  SS.amplifiedRegion(i,find(pos_f(sameSeq{i}(1),:))) = find(pos_f(sameSeq{i}(1),:));
  
  if ~isempty(find(curr_seq-int2nt(unpack_seqs(SS.Sequence_packed_matrix(i,:),readLength*2*6,64))))
    disp('problem')
    pause
  end
  
    
  origSeq_tmp = zeros(length(sameSeq{i}),readLength*6*2);
  SS.Header_uni{i} = cell(1,length(sameSeq{i}));
  for j=1:length(sameSeq{i})
    
    % check the group is indeed the same over the relevant regions
    aa_f = find(pos_f(sameSeq{i}(j),:));curr_ind_f = pos_f(sameSeq{i}(j),aa_f);
    aa_r = find(pos_r(sameSeq{i}(j),:));curr_ind_r = pos_r(sameSeq{i}(j),aa_r);
    for k=1:length(aa_f)
      origSeq_tmp(j,(aa_f(k)-1)*readLength+1:aa_f(k)*readLength) = S.Sequence_uni{sameSeq{i}(j)}(curr_ind_f(k):curr_ind_f(k)+readLength-1);
      origSeq_tmp(j,basicOffSet+(aa_f(k)-1)*readLength+1:basicOffSet+aa_f(k)*readLength) = S.Sequence_uni{sameSeq{i}(j)}(curr_ind_r(k):curr_ind_r(k)+readLength-1);
    end
    
    SS.Header_uni{i}{j} = S.Header_uni{sameSeq{i}(j)};
  end
  
  if size(unique(origSeq_tmp,'rows'),1)>1
    disp('problem')
    pause
  end
  
  if isempty(SS.amplifiedRegion(i,:))
    disp('problem')
    pause
  end
  
  % check all are amplified in the same regions
  curr_amp = pos_f(sameSeq{i},:);
  curr_amp(find(curr_amp)) = 1;
  if size(unique(curr_amp,'rows'),1)>1
    disp('problem')
    pause
  end
  %%%%%%%%%
  
end

SS.Sequence_packed = cell(1,numBact);
for i=1:numBact
  if mod(i,1000)==0
    i
  end
  SS.Sequence_packed{i} = SS.Sequence_packed_matrix(i,:);
end

%SS = rmfield(SS,'Sequence_packed_matrix');


save([bcs_data_path,'full16S/data_primersAmit_2MM_databaseOf2012_readLength_76_noACGT_onlySS'],'SS') 


% matrix 

%
%pause
clear values M 
for reg=1:6
  reg
  % output the bacteria in each region
  [PrimerReadBySpeciesMat_tmp Primer_kmers_packed_tmp bactInRegion_tmp]=oneRegionDataForConcatenatedData(reg,readLength,SS);
  
  [junk,i1] = intersect(1:numBact,bactInRegion_tmp);
  disp('is this true?')
  M{reg} = sparse(size(PrimerReadBySpeciesMat_tmp,1),numBact);
  M{reg}(:,i1) = PrimerReadBySpeciesMat_tmp(:,bactInRegion_tmp);
  values{reg} = Primer_kmers_packed_tmp;
end

save([bcs_data_path,'full16S/data_primersAmit_2MM_databaseOf2012_readLength_76_noACGT'],'sameSeq','seq_readLength76','M','values') 

%whos  sameSeq   seq_readLength76   M   values 

%check - overlap bwtwqeen region
for i=1:5
  for j=i+1:6
    [junk i1 i2]=intersect(values{i},values{j},'rows');
    if ~isempty(junk)
      disp('overlap!!!')
      pause
    end
  end
end

