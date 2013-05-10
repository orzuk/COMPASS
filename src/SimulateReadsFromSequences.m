% Simulate reads from a set of sequences given as
%
% Input:
% readLength - length of reads to extract (currently only single-end reads)
% packed_seqs - array of input sequences (in packed form)
% seqs_lens - length of input sequences (in nucleotides)
% num_reads - total number of reads to simulate
% seqs_weights_vec - fraction of each sequence in simulated reads
% read_error_substitution_table - 4x4xreadLength error substitution matrix on extracted reads (optional)
%
% Output:
% noisy_kmers_packed - simulated kmers (with noise model) 
% clean_kmers_packed - simulated clean kmers (before noise is added) 
% kmer_inds - index of each kmer in the species list 
%
function [noisy_kmers_packed clean_kmers_packed kmer_inds] = ...
    SimulateReadsFromSequences(readLength, packed_seqs, seqs_lens, num_reads, ...
    seqs_weights_vec, read_error_substitution_table)

num_seqs = size(packed_seqs, 1);
if(length(seqs_lens) == 1)
    seqs_lens = repmat(seqs_lens, num_seqs, 1);
end

% First randomize positions (internal function)
[kmer_seq_inds, kmer_positions_inds] = ...
    RandomizeReadsPoisitions(seqs_weights_vec, seqs_lens, num_reads);

% Next extract reads according to positions
[clean_kmers_packed kmer_inds] = ...
    extract_sub_kmers(packed_seqs, seqs_lens, readLength, 0, 0, ...
    kmer_seq_inds, kmer_positions_inds);

% Finally add noise according to noise model
if(exist('read_error_substitution_table', 'var'))
    noisy_kmers_packed = add_noise_to_kmers(readLength, num_reads, ...
        clean_kmers_packed, reshape(read_error_substitution_table, readLength, 4*4));
else
    noisy_kmers_packed = clean_kmers_packed;
end

% Internal function: randomize sequences and positions according to sequence weights
% Currently supported only uniform model along the sequence
function [kmer_seq_inds, kmer_positions_inds] = ...
    RandomizeReadsPoisitions(seqs_weights_vec, seqs_lens, num_reads)

kmer_seq_inds = vec2column(weighted_rand(vec2row(seqs_weights_vec), num_reads)); % randomize sequence ids
kmer_positions_inds = ceil(rand(num_reads,1) .* seqs_lens(kmer_seq_inds));








