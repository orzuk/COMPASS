% Prepare the mixing matrix from a set of sequences and read length
%
% Input:
% readLength - length of reads
% packed_seqs - sequences array (should be in packed form. 2bits per nucleotide)
% seqs_lens - length in nucleotide of each sequence
% matlab_word_size - % default is 64
% kmer_seq_inds - NEW! indices of sequences to use
% kmer_positions_inds - NEW! indices of positions to use
%
% Output:
% ReadBySpeciesMat - matrix with 1 if read i appears in sequence j (in sparse format)
% kmers_packed - read sequences for each row of A 
%
function [ReadBySpeciesMat kmers_packed] = BuildMixingMatrixFromSequences( ...
    readLength, packed_seqs, seqs_lens, matlab_word_size, kmer_seq_inds, kmer_positions_inds)

sprintf('Prepare kmers')
AssignGeneralConstants; % set constants including word size 
if(~exist('matlab_word_size', 'var'))
    matlab_word_size = 64; % set default: unix 64
end
if(matlab_word_size == 32)
    machine_word_str = 'uint32'
else
    machine_word_str = 'uint64'
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
hash_flag = 0; % don't use hasing (not working) 
debug_prints = 0; % print debug statements in C function 

sprintf('Extract kmers') 
whos packed_seqs 
if(exist('kmer_seq_inds', 'var') && (~isempty(kmer_seq_inds)))
    [kmers_packed kmers_inds] = ...
        extract_sub_kmers(packed_seqs, seqs_lens, readLength, unique_flag, [hash_flag debug_prints], ...
        kmer_seq_inds, kmer_positions_inds); % Call c-function (no hash!!!) with pre-specified locations
else                     
    [kmers_packed kmers_inds] = ...
        extract_sub_kmers(packed_seqs, seqs_lens, readLength, unique_flag, [hash_flag debug_prints]); % Call c-function (no hash!!!)
end
ReadBySpeciesMat = sparse(kmers_inds(:,1), kmers_inds(:,2), 1); % generate a sparse reads matrix



%keyboard
%return;

% % % run_noam_version = 0;
% % % if(run_noam_version) % Below is Noam's version
% % %     o = 1:num_species;
% % %     for i=1:num_species % loop on all species
% % %         if mod(i,1000)==0
% % %             i
% % %         end
% % %
% % %         seq{i} = char(55*ones(length(packed_seqs{o(i)})-readLength,readLength)); % what's these 3 lines dowing?
% % %         for j=0:length(packed_seqs{o(i)})-readLength
% % %             seq{i}(j+1,:) = packed_seqs{o(i)}(j+1:j+readLength);
% % %         end
% % %     end
% % %     clear packed_seqs  % get rid of input to save space
% % %
% % %
% % %     for i=1:num_spe
% % %
% % %         ln = zeros(1,length(o));
% % %         for i=1:length(o)
% % %             ln(i) = size(seq{i},1);
% % %         end
% % %         len = sum(ln);
% % %
% % %         k = 0;
% % %         seqO = char(55*ones(len,readLength));
% % %         posO = zeros(len,1);
% % %         for i=1:length(o)
% % %             seqO(k+1:k+length(seq{i}),:) = seq{i};
% % %             posO(k+1:k+length(seq{i})) = i*ones(length(seq{i}),1);
% % %             k = k+length(seq{i});
% % %         end
% % %         clear seq
% % %
% % %         %%%%%%%
% % %         % find the reads in each of these
% % %
% % %         % output the data
% % %         [values,currInds,positions,leng] = myUnique(seqO,posO);
% % %         % values = the unique rows in Esq
% % %         % currInds = the indices of the unique rows in seqO
% % %         % positions = the BAC which holds the sequence and the number of times it does.
% % %         % leng = the total number of BAC which hold the sequence
% % %
% % %         clear seqO currInds posO
% % %
% % %         [vals_leng inds_leng] = get_duplicates(leng);
% % %         % inds_leng = indices of unique values of leng - the interesting one is leng==1
% % %
% % %         clear leng
% % %
% % %         %pause
% % %         keep = 1:size(values,1);
% % %
% % %         lenb = 0;
% % %         for i=1:length(inds_leng)
% % %             i
% % %             i1 = 1:length(inds_leng{i});
% % %             [junk,i11,i22] = intersect(inds_leng{i}(i1),keep);
% % %             b = i22'*ones(1,vals_leng(i));
% % %             lenb = lenb+length(b(:));
% % %         end
% % %         % if works - this can be changed
% % %
% % %         lenKeep = length(keep);
% % %         put = zeros(lenb,3);
% % %         k = 0;
% % %         for i=1:length(inds_leng)
% % %             i
% % %             i1 = 1:length(inds_leng{i});
% % %             [junk,i11,i22] = intersect(inds_leng{i}(i1),keep);
% % %             if length(i11)~=length(i1)
% % %                 disp('error')
% % %                 pause
% % %             end
% % %             b = i22'*ones(1,vals_leng(i));
% % %
% % %             dat = cell2mat(positions(inds_leng{i}));
% % %             curr_pos = dat(1:2:end,:);curr_pos = curr_pos(i1,:);
% % %             curr_num = dat(2:2:end,:);curr_num = curr_num(i1,:);
% % %
% % %             if ~isempty(find(curr_pos==0))
% % %                 i
% % %                 pause
% % %             end
% % %             curr_dat = [b(:) curr_pos(:) curr_num(:)];
% % %             put(k+1:k+size(curr_dat,1),:) = curr_dat;
% % %             k = k+size(curr_dat,1);
% % %
% % %         end
% % %
% % %         clear positions inds_leng d
% % %         % change this - to include values themselves and the work on the recuded number of rows
% % %
% % %
% % %         part = 1:10^5:lenb; % divide to blocks
% % %         if part(end)<lenb
% % %             part(end+1) = lenb;
% % %         end
% % %         part = [part,lenb+1];
% % %
% % %         ReadBySpeciesMat = spalloc(length(keep),length(o),lenb);
% % %         for i=1:length(part)-1
% % %             i
% % %             inWhichBAC = spalloc(length(keep),length(o),lenb);
% % %             do = part(i):part(i+1)-1;
% % %             inWhichBAC((put(do,2)-1)*lenKeep+put(do,1)) = put(do,3);
% % %             ReadBySpeciesMat = ReadBySpeciesMat+inWhichBAC;
% % %         end
% % %         clear put
% % %         nnz(ReadBySpeciesMat)
% % %         clear inWhichBAC keep
% % %
% % %         for i=1:length(o) % normalize columns: each bacteria is normalized by its length
% % %             ReadBySpeciesMat(:,i) = ReadBySpeciesMat(:,i)/ln(i);
% % %         end
% % %         %keyboard
% % %     end %
% % % end % if run noam's version




% Prepare the mixing matrix from a set of sequences and read length. Here
% input also subset of positions
%
% Input:
% readLength - length of reads
% packed_seqs - sequences array (should be in packed form. 2bits per nucleotide)
% seqs_lens - length in nucleotide of each sequence
%
% Output:
% ReadBySpeciesMat - matrix with 1 if read i appears in sequence j (in sparse format)
% kmers_packed - read sequences for each row of A (isn't it the input??)
%
% % % function [ReadBySpeciesMat kmers_packed] = BuildMixingMatrixFromSequences( ...
% % %     readLength, packed_seqs, seqs_lens, matlab_word_size)

