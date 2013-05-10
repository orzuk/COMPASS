disp('repeat this with 2MM and the whole database')

clear

bcs_data_path = '~/CS/BAC/'

primers_locations_file = fullfile(bcs_data_path, 'AmitPrimers/primers_position_few_regions.mat');

% database:
seq_database_file = fullfile(bcs_data_path, '/full16S/bac16s_full_without_ambiguous'); % Database with all sequences
len_uni = load([bcs_data_path,'datNoNonACGT/keyNoNonACGT'])

P = load(primers_locations_file);   % load primer locations

load([bcs_data_path,'/data_before_right_edge_correction/listOfDoubles']);
list = 1:455055;
list(d) = [];
primer_fields = fields(P);
for i=1:length(primer_fields)
    eval_str = ['P_no_nonACGT.',primer_fields{i} ' = ','P.',primer_fields{i},'(list,:);'];    
    eval(eval_str); 
end
save([bcs_data_path,'full16S/primerDataAmitPrimers_full_without_ambiguous'],'P_no_nonACGT')




if 1==2 % done only once
  S = load(seq_database_file); 
  clear SS
  SS.len_uni = S.len_uni(inds);
  SS.Header_uni = S.Header_uni(inds);
  for i=1:length(inds)
    SS.Sequence_packed{i} = pack_seqs(S.Sequence_uni{inds(i)},64);
  end
  save ~/CS/BAC/full16S/Sequence_packed_full_without_ambiguous SS
  
end


% check databaes - with amb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read length 26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
bcs_data_path = '~/CS/BAC/'
seq_database_file = fullfile(bcs_data_path, '/full16S/bac16s_full_without_ambiguous'); % Database with all sequences
load([bcs_data_path,'full16S/primerDataAmitPrimers_full_without_ambiguous'])

readLength = 26; % set read length  - this corresponds to 100 !!!!
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
  S.len_uni = vec2column(len_uni.len_uni); 
  
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
  
  [tmp_seq_readLength26,origInd] = sortrows([seq_f,seq_r]);


  [u_f,i1_f] = unique(tmp_seq_readLength26,'rows','first');
  [u_l,i1_l] = unique(tmp_seq_readLength26,'rows','last');
  
  seq_readLength26 = u_f;
  
  % find the bacteria that are part of this 
  sameSeq = cell(length(i1_f),1);
  group = zeros(length(S.Sequence_uni),1);
  for i=1:length(i1_f)
    sameSeq{i} = origInd(i1_f(i):i1_l(i));
    group(sameSeq{i}) = i;
  end
  
  
  save([bcs_data_path,'full16S/data_primersAmit_0MM_readLength_26'],'sameSeq','seq_readLength26','group') 
end % create the reads 





% matrix 

% prepare matrix for each region
load  ~/CS/BAC/full16S/Sequence_packed_full_without_ambiguous

% check whether each region is unique

clear M value
inds = origInd(i1_f); % take a representative
for reg=1:6
  % output the bacteria in each region
  [PrimerReadBySpeciesMat_tmp Primer_kmers_packed_tmp bactInRegion_tmp]=oneRegionData(reg,inds,seq_ind,pos_f,pos_r,readLength,SS);
  
  [junk,i1] = intersect(inds,bactInRegion_tmp);
  
  M{reg} = sparse(size(PrimerReadBySpeciesMat_tmp,1),length(inds));
  M{reg}(:,i1) = PrimerReadBySpeciesMat_tmp(:,bactInRegion_tmp);
  values{reg} = Primer_kmers_packed_tmp;
end

save([bcs_data_path,'full16S/data_primersAmit_0MM_readLength_26'],'sameSeq','seq_readLength26','group','M','values') 

%check - overlap bwtwqeen region
for i=1:5
  for j=i+1:6
    junk=intersect(values{i},values{j},'rows');
    if ~isempty(junk)
      disp('overlap!!!')
    end
  end
end


