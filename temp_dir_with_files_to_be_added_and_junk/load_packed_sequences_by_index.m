% Load packed sequences of different species from set of files
%
% Input:
% basicSeqKey - key mapping sequences indices to correct files
% seqs_inds - vector of indices of each sequence we want to load
%
% Output
% packed_seqs - cell array of sequences in packed form
% seqs_lens - length of each sequences (in nucleotides)
%
function [packed_seqs seqs_lens] = load_packed_sequences_by_index(basicSeqKey, seqs_inds)

load(basicSeqKey, 'positionInPart','len_uni'); % load key mapping sequences to files
basicSeqNameDir = dir_from_file_name(basicSeqKey); % directory with all sequences stored in separate files
num_seqs = length(seqs_inds);
seqs_lens = len_uni(seqs_inds); % get lengths

packed_seqs = cell(num_seqs,1);
for i=1:num_seqs % loop on sequences and load relevant files
    clear seq_*
    load([basicSeqNameDir, 'seq_part_',num2str(positionInPart(seqs_inds(i)))], ...
        ['seq_',num2str(seqs_inds(i))]);
    w = ['packed_seqs{i} = seq_',num2str(seqs_inds(i)),'{1};'];
    eval(w);
    
end
clear seq_* positionInPart
