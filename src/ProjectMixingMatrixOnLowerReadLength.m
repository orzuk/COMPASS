% Project a large mixing matrix with large L on a smaller L. 
% Problem: the function is working but is slower than re-computing the
% matrix from the sequences - how to make this faster? 
% 
% Input:
% original_readLength - length of reads
% projected_readLength - length of new (shorter) reads 
% packed_seqs - sequences array (should be in packed form. 2bits per nucleotide)
% seqs_lens - length in nucleotide of each sequence
% matlab_word_size - % default is 64
%
% Output:
% ProjectedReadBySpeciesMat - matrix with 1 if read i appears in sequence j using the new read length L (in sparse format)
% projected_kmers_packed - read sequences for each row of A
%
function [ProjectedReadBySpeciesMat projected_kmers_packed] = ProjectMixingMatrixOnLowerReadLength( ...
    ReadBySpeciesMat, kmers_packed, original_readLength, projected_readLength, packed_seqs, seqs_lens, matlab_word_size)

if(~exist('matlab_word_size', 'var'))
    matlab_word_size = 64; % set default: unix 64
end
if(matlab_word_size == 32)
    machine_word_str = 'uint32';
else
    machine_word_str = 'uint64';
end
% if(iscell(packed_seqs)) % convert to an array (this is OK if lengths are not TOO different)
%     num_species = length(packed_seqs);
%     seq_lens_in_words = length_cell(packed_seqs);
%     packed_seqs_mat = zeros(num_species, max(seq_lens_in_words), machine_word_str);
%     for i=1:num_species
%         packed_seqs_mat(i,1:seq_lens_in_words(i)) = packed_seqs{i};
%     end
%     packed_seqs = packed_seqs_mat; clear packed_seqs_mat; % save space
% end

% % % % % % kmers_unpacked = unpack_seqs(kmers_packed, original_readLength); % This part is too heavy
% % % % % % kmers_unpacked = kmers_unpacked(:, 1:projected_readLength); % take only prefix with shorter read length 
% % % % % % [projected_kmers_packed I J] = unique(kmers_unpacked, 'rows'); % First just get the unique kmers - all the prefixes of the original kmers 
% % % % % % projected_kmers_packed = pack_seqs(projected_kmers_packed); % pack back 
num_kmers = size(kmers_packed, 1);

[projected_kmers_packed kmers_inds2] = ... % This one causes matlab to crash
    extract_sub_kmers(kmers_packed, repmat(original_readLength, num_kmers, 1), projected_readLength, 0, 0, ...
    1:num_kmers, ones(num_kmers, 1)); % Call c-function (no hash!!!) with pre-specified locations
[projected_kmers_packed I J] = unique(projected_kmers_packed, 'rows'); % First just get the unique kmers - all the prefixes of the original kmers 

[kmer_inds seqs_inds num_kmers_in_seqs] = find(ReadBySpeciesMat);  % get indices and values of non-zero entries in sparse matrix
kmer_inds = J(kmer_inds); 

%ProjectedReadBySpeciesMat = sparse(kmer_inds, seqs_inds, num_kmers_in_seqs); % generate a sparse reads matrix
%return; 

num_seqs = length(packed_seqs);  % Correct for tail 
for i=1:num_seqs
    unpacked_seqs{i} = unpack_seqs(packed_seqs{i}, seqs_lens(i)); % correct for tail 
    unpacked_seqs{i} = unpacked_seqs{i}((seqs_lens(i)-original_readLength+2):end); % assume all are long enough
end
packed_seqs = pack_seqs(unpacked_seqs); 
[TailReadBySpeciesMat tail_kmers_packed] = BuildMixingMatrixFromSequences(projected_readLength, packed_seqs, ...
    repmat(original_readLength-1, num_seqs, 1), matlab_word_size); 

[tail_kmer_inds tail_seqs_inds tail_num_kmers_in_seqs] = find(TailReadBySpeciesMat);  % get indices and values of non-zero entries in sparse matrix


[union_kmers_packed] = union(projected_kmers_packed, tail_kmers_packed, 'rows'); % Unite two matrices 
[~, ~, I_union] = intersect(projected_kmers_packed, union_kmers_packed, 'rows'); 
[~, ~, J_union] = intersect(tail_kmers_packed, union_kmers_packed, 'rows'); 

union_kmer_inds = [vec2row(I_union(kmer_inds)) vec2row(J_union(tail_kmer_inds))]'; 
union_seqs_inds = [seqs_inds' tail_seqs_inds']'; 
union_num_kmers_in_seqs = [num_kmers_in_seqs' tail_num_kmers_in_seqs']'; 

ProjectedReadBySpeciesMat = sparse(union_kmer_inds, union_seqs_inds, union_num_kmers_in_seqs); % generate a sparse reads matrix
projected_kmers_packed = union_kmers_packed;








%projected_kmer_inds = projected_kmers_packed(kmer_inds); % new indices with projected kmers 
%[new_inds new_I new_J] = unique([projected_kmer_inds seqs_inds], 'rows');


%accum_num_kmers = accumarray(new_J, num_kmers_in_seqs);



% Extract 



