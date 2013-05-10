function [PrimerReadBySpeciesMat Primer_kmers_packed bactInRegion]=oneRegionData(reg,inds,seq_ind,pos_f,pos_r,readLength,SS)

%keyboard
p_seq_inds = seq_ind(inds,reg);
p_position_inds_f = pos_f(inds,reg);
p_position_inds_r = pos_r(inds,reg);
a = find(p_seq_inds==0);
p_seq_inds(a) = [];
p_position_inds_f(a) = [];
p_position_inds_r(a) = [];

p_seq_inds = vec2row(p_seq_inds(:));
p_position_inds_f = vec2row(p_position_inds_f(:));
p_position_inds_r = vec2row(p_position_inds_r(:));

bactInRegion = unique(p_seq_inds);

p_seq_inds = [p_seq_inds, p_seq_inds];
p_position_inds = [p_position_inds_f,p_position_inds_r];


%keyboard
[PrimerReadBySpeciesMat Primer_kmers_packed] = BuildMixingMatrixFromSequencesRegions( ...
    readLength, SS.Sequence_packed, SS.len_uni, [], p_seq_inds,p_position_inds); % matrix using only pre-specified locations
