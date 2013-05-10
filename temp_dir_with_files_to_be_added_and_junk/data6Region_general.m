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

