% Added new option: build reads-by-species matrix where only certain positions along the sequences are valid
% function test_new_build_matrix
AssignGeneralConstants;

test_build_matrix=0;
test_projection=1;

switch machine
    case PC
        bcs_data_path = 'C:\research\compressed_sensing\metagenomics\nextgen\data\'
        bcs_path = '../../compressed_sensing/metagenomics/next_gen';
    case UNIX
        bcs_data_path = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/data';
end
primers_locations_file = fullfile(bcs_data_path, 'AmitPrimers', 'primers_position_few_regions.mat');
seq_database_file = fullfile(bcs_data_path, 'bacteria_s16_data_uni.mat'); % Database with all sequences
primers_matrix_output_file = fullfile(bcs_data_path, 'AmitPrimers', 'PrimersReadsBySequencesMatrix.mat');
readLength = 50; % set read length


% fullfile(bcs_path, 'data'); % here data is stored
% bcs_figs_path = fullfile(bcs_data_path, 'figs'); % here we save output figures
% bcs_output_path = fullfile(bcs_path, 'results'); % here put simulations results
% reads_data_file = fullfile(bcs_data_path, 'bac_sim_100outof10000_10e6_reads.mat'); % 'bac_sim_seqdata_100outof10000_10e6_reads_similarbac_n0p5_freq_3.mat'); % file with all reads
% comparison_output_file = '~/CS/BAC/res_10e6reads_100outof450000_similar_usingIteration'; % where to save comparison



P = load(primers_locations_file);   % load primer locations
S = load(seq_database_file); S.Sequence_uni = vec2column(S.Sequence_uni); S.len_uni = vec2column(S.len_uni);
max_ind = length(S.len_uni); % .100; % for debug
S = struct_by_inds(S, 1:max_ind); % for debug: save time run only a few sequences


primer_fields = fields(P); primer_seq_inds = []; primer_positions_inds = [];
for i=1:length(primer_fields)
    eval_str = ['primer_seq_inds = [primer_seq_inds 1:length(P.' primer_fields{i} ')];'];
    eval(eval_str); eval(eval_str); % twice (once for forward and once for reverse)
    eval_str = ['primer_positions_inds = [primer_positions_inds vec2row(P.' primer_fields{i} '(:,1))];'];  % forward strand
    eval(eval_str);
    eval_str = ['primer_positions_inds = [primer_positions_inds vec2row(P.' primer_fields{i} '(:,2)) - readLength];'];  % reverse strand
    eval(eval_str);
end

good_inds = find(primer_positions_inds > 0); % get rid of zero indices
good_inds = intersect(good_inds, find(primer_seq_inds <= max_ind)); % temp for debug
primer_positions_inds = primer_positions_inds(good_inds);
primer_seq_inds = primer_seq_inds(good_inds);


if(test_build_matrix)
    if(max_ind <= 1000) % make full matrix for comparison
        [ReadBySpeciesMat kmers_packed] = BuildMixingMatrixFromSequences( ...
            readLength, S.Sequence_packed, S.len_uni); % full matrix
    end
    
    %
    [PrimerReadBySpeciesMat Primer_kmers_packed] = BuildMixingMatrixFromSequences( ...
        readLength, S.Sequence_packed, S.len_uni, [], primer_seq_inds, primer_positions_inds); % matrix using only pre-specified locations
    
    save(primers_matrix_output_file , 'PrimerReadBySpeciesMat', 'Primer_kmers_packed'); % save output to file
    
    small_seqs = {'AAAAAAAAGAGGGC', 'CCCCCCCTTTTTCGT'};
    [small_packed_seqs small_packed_seqs_lens] = pack_seqs(small_seqs, 32);
    
    [SmallReadBySpeciesMat small_kmers_packed] = BuildMixingMatrixFromSequences( ...
        10, small_packed_seqs, small_packed_seqs_lens); % matrix using only pre-specified locations
    full(SmallReadBySpeciesMat)
    int2nt(unpack_seqs(small_kmers_packed, 10))
end % if test_build_matrix




if (test_projection) % here test projection
    small_readLength = 10;
    [PrimerReadBySpeciesMat Primer_kmers_packed] = BuildMixingMatrixFromSequences( ...
        readLength, S.Sequence_packed, S.len_uni, [], primer_seq_inds, primer_positions_inds); % matrix using only pre-specified locations
    [SmallPrimerReadBySpeciesMat SmallPrimer_kmers_packed] = BuildMixingMatrixFromSequences( ...
        small_readLength, S.Sequence_packed, S.len_uni, [], primer_seq_inds, primer_positions_inds); % matrix using only pre-specified locations
    
    [ProjectedSmallPrimerReadBySpeciesMat ProjectedSmallPrimer_kmers_packed] = ProjectMixingMatrixOnLowerReadLength( ...
        PrimerReadBySpeciesMat, Primer_kmers_packed, readLength, small_readLength);  % project matrix
    
end % test projection


