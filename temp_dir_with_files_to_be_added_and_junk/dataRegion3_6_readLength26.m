%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read length 26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
bcs_data_path = '~/CS/BAC/';
seq_database_file = fullfile(bcs_data_path, '/full16S/bac16s_full_without_ambiguous'); % Database with all sequences
load([bcs_data_path,'full16S/primerDataAmitPrimers_full_without_ambiguous'])

readLength = 26; % set read length  - this corresponds to 100 !!!!
AmitBarcodeOffset = 18; % start read from this.


primer_fields = fields(P_no_nonACGT); 
seq_ind = zeros(size(P_no_nonACGT.pos_pr1,1),length(primer_fields));
pos_f = seq_ind;
pos_r = seq_ind;
for i=[3:6]%1:length(primer_fields)
  
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
  
  
  
end % create the reads 





% matrix 

% prepare matrix for each region
load  ~/CS/BAC/full16S/Sequence_packed_full_without_ambiguous

% check whether each region is unique

clear M value
[group_val,inds]= unique(group);  % take a representative

1
pause
for reg=3:6
  % output the bacteria in each region
  [PrimerReadBySpeciesMat_tmp Primer_kmers_packed_tmp bactInRegion_tmp]=oneRegionData(reg,inds,seq_ind,pos_f,pos_r,readLength,SS);
  
  [junk,i1] = intersect(inds,bactInRegion_tmp);
  
  M{reg} = sparse(size(PrimerReadBySpeciesMat_tmp,1),length(inds));
  M{reg}(:,i1) = PrimerReadBySpeciesMat_tmp(:,bactInRegion_tmp);
  values{reg} = Primer_kmers_packed_tmp;
end

save([bcs_data_path,'full16S/data_primersAmit_0MM_readLength_26_noReg12'],'sameSeq','seq_readLength26','group','M','values','inds') 

%check - overlap bwtwqeen region
for i=1:5
  for j=i+1:6
    junk=intersect(values{i},values{j},'rows');
    if ~isempty(junk)
      disp('overlap!!!')
    end
  end
end


