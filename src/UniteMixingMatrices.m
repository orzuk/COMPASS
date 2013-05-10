% Unite several mixing matrices obtained from different sets of species 
%
% Input:
% ReadBySpeciesMat - a cell-array of (sparse) matrices
% kmers_packed - a cell-array of kmers vectors (packed form) 
% readLength - length of reads
% packed_seqs - sequences array (should be in packed form. 2bits per nucleotide)
% seqs_lens - length in nucleotide of each sequence
% matlab_word_size - % default is 64
% kmer_seq_inds - NEW! indices of sequences to use
% kmer_positions_inds - NEW! indices of positions to use
%
% Output:
% UnitedReadBySpeciesMat - matrix with 1 if read i appears in sequence j for all species (in sparse format)
% united_kmers_packed - read sequences for each row of A 
%
function [UnitedReadBySpeciesMat united_kmers_packed] = UniteMixingMatrices( ...
    ReadBySpeciesMat, kmers_packed, ... % vector of the original matrices to unite 
    readLength, packed_seqs, seqs_lens, matlab_word_size, kmer_seq_inds, kmer_positions_inds)

AssignGeneralConstants; % set constants including word size 
if(~exist('matlab_word_size', 'var'))
    matlab_word_size = 64; % set default: unix 64
end
if(matlab_word_size == 32)
    machine_word_str = 'uint32';
else
    machine_word_str = 'uint64';
end


num_matrices = length(ReadBySpeciesMat); % get number of matrices;
UnitedReadBySpeciesMat = ReadBySpeciesMat{1};
united_kmers_packed = kmers_packed{1};



for i=2:num_matrices % perform union iteratively
    united_kmers_packe = union(united_kmers_packed, kmers_packed{i}, 'rows');  % unite the kmers 
end
intersect_inds = cell(num_matrices, 1); 
for i=1:num_matrices % get indices of each small matrix in the big matrix 
    [~, intersect_inds{i}] = intersect(united_kmers_packed, kmers_packed{i}); 
end


if(iscell(packed_seqs)) % convert to an array (this is OK if lengths are not TOO different)
    num_species = length(packed_seqs);
    seq_lens_in_words = length_cell(packed_seqs);
    packed_seqs_mat = zeros(num_species, max(seq_lens_in_words), machine_word_str);
    for i=1:num_species
        packed_seqs_mat(i,1:seq_lens_in_words(i)) = packed_seqs{i};
    end
    packed_seqs = packed_seqs_mat; clear packed_seqs_mat; % save space
end

unique_flag = 1; % flag saying to return only unique kmers

UnitedReadBySpeciesMat = sparse(kmers_inds(:,1), kmers_inds(:,2), 1); % generate a sparse reads matrix



