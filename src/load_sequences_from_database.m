% Load sequences of specified bacteria from database
%
% Input:
% database_type - how is database stored (default: in blocks)
% database_dir - database directory
% database_keyfile - database keyfile (or database file itself)
% seq_inds - indices of species to extract
%
% Output:
% packed_seqs - sequences of extracted species (packed form)
% seqs_lens - lengths of extracted sequences (in nucleotides)
%
function [packed_seqs seqs_lens] = load_sequences_from_database(database_type, database_dir, database_keyfile, seq_inds)


if(~exist('database_type', 'var') || isempty(database_type))
    database_type = 'partitioned_to_blocks'; % blocks (Noam's format)
end

switch database_type
    case {'partitioned_to_blocks', 'blocks'} % Database partitioned to blocks from Noam (fast)
        
        % load sequences of bacteria in the mixture
        load(database_keyfile, 'positionInPart', 'len_uni');
        seqs_lens = len_uni(seq_inds); clear len_uni
        num_seqs = length(seq_inds);
        packed_seqs = cell(num_seqs,1);
        for i=1:num_seqs
            clear seq_*
            load(fullfile(database_dir, 'seq_part_',num2str(positionInPart(seq_inds(i)))), ...
                ['seq_',num2str(seq_inds(i))]);
            eval_str = ['packed_seqs{i} = seq_',num2str(seq_inds(i)),'{1};'];
            eval(eval_str);
        end
        clear seq_* positionInPart
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case {'one-big-file'} % My version (slow since we load all sequences) 
        load(database_keyfile, 'packed_seqs', 'seqs_lens');
        packed_seqs = packed_seqs(seq_inds,:);
        seqs_lens = seqs_lens(seq_inds);
end






