% Iterate over ...
%
% Input:
% uniqueReads - set of reads (unique)
% uniqueReads_counts - counts of reads (how many copies of each read)
% species_packed_seqs - what's that? the sequences of all bacterial species
% species_names - names of bacterias
% readLength - length of reads
% block_size - size of blocks (solve each block seperately)
%
% Output:
% X - vector of all solutions (frequency of each bacteria) - concatenation of block solutions
%
function X = iterate(uniqueReads, uniqueReads_counts, ...
    species_packed_seqs, species_names, readLength, ...
    algorithm_str, block_size, algorithm_parameters)

if(~exist('block_size', 'var') || isempty(block_size)) % set default block size 
    block_size = 1000;
end
num_species = size(species_names,1); % n is ???
X = zeros(1,num_species); % Solution vector

switch algorithm_str % choose which algorithm to run
    case 'block_by_block' % solve on blocks using cvs
        %        block_size = algorithm_parameters;
        [block_starts, block_ends, ~, num_blocks] = ...
            divide_region_to_blocks(1, num_species, block_size); % divide to blocks
        rand_P = randperm(num_species);
        firstFlag = 1;
        x_blocks = cell(num_blocks, 1);
        for i=1:num_blocks % solve seperately for each block of 1000
            block_seqs = species_packed_seqs(rand_P(block_starts(i):block_ends(i)));   % create a part of the matrix
            [ReadBySpeciesMat values] = BuildMixingMatrixFromSequences(readLength,block_seqs);
            [x_blocks{i},firstFlag] = SolveMixingMatrixFromReads(firstFlag, ReadBySpeciesMat, values, ...
                uniqueReads,uniqueReads_counts,0);            % solve for one block (1000)
        end
        
        for i=1:num_blocks % Concatenate results: this is not normalized
            if ~isempty(x_blocks{i})
                X(rand_P(block_starts(i):block_ends(i))) = x_blocks{i};
            end
        end
    case 'stochastic_gradient'
        num_iters = algorithm_parameters;
        X = l2cd(ReadBySpeciesMat, b, num_iters); 
    case 'gradient_decent'
        num_iters = algorithm_parameters;
        X = l2mu(ReadBySpeciesMat, b, num_iters);         
end % switch which algorithm to use
