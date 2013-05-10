%clear
AssignGeneralConstants;
readLength = 50; % read length
num_reads_used = 10^6; % don't use all reads (make everything faster)
num_species_in_database = 1000; % use a subset of bacteria (make everything faster)
num_species_in_mixture = 100;  % how complicated is the mixture

switch machine
    case PC
        plot_summary_statistics = 1; % plot statistics for database
    case UNIX
        plot_summary_statistics = 0; % plot statistics for database
end
block_size = 1000; % size of blocks to use when solving
algorithm_str = 'block_by_block'; % how to solve 
algorithm_str_vec = {'block_by_block'};

cost_function_str = 'L2'; % what are we trying to optimize 
switch machine
    case PC
        bcs_path = '../../compressed_sensing/metagenomics/next_gen';
        %        output_data_path = '~/CS/BAC/iterations100outof450000_based10e6_similar';
    case UNIX
        bcs_path = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
        %        output_data_path = '~/CS/BAC/iterations100outof450000_based10e6_similar';
end
bcs_data_path = fullfile(bcs_path, 'data'); % here data is stored
bcs_figs_path = fullfile(bcs_data_path, 'figs'); % here we save output figures
bcs_output_path = fullfile(bcs_path, 'results'); % here put simulations results
seq_database_file = fullfile(bcs_data_path, 'bacteria_s16_data_uni.mat'); % Database with all sequences
reads_data_file = fullfile(bcs_data_path, 'bac_sim_100outof10000_10e6_reads.mat'); % 'bac_sim_seqdata_100outof10000_10e6_reads_similarbac_n0p5_freq_3.mat'); % file with all reads
comparison_output_file = '~/CS/BAC/res_10e6reads_100outof450000_similar_usingIteration'; % where to save comparison


t_read_reads = cputime;
load(reads_data_file, 'reads_data_packed','w','Header_mix'); % load reads , weights (answer), header_mix (bacterial names)
if(num_reads_used < length(reads_data_packed))
    P = randperm(length(reads_data_packed));
    reads_data_packed = reads_data_packed(P(1:num_reads_used),:); % take a subset of reads
    save(fullfile(bcs_data_path, 'bac_sim_100outof10000_10e6_reads.mat'), ...
        'reads_data_packed','w','Header_mix');
end
% [uniqueReads,uniqueReads_inds] = myUnique(reads_data); % perform unique on all reads
[uniqueReads,uniqueReads_counts] = unique_with_counts(reads_data_packed, 'rows');  % perform unique on all reads
clear reads_data_packed; % save memory
t_read_reads = cputime - t_read_reads

t_read_database = cputime;
if(~exist('Sequence_packed', 'var'))
    load(seq_database_file, ...
        'Header_uni', 'len_uni', 'Sequence_packed');  % load database (in packed form). This is still too long (~1/2 minute). Need to save database in blocks (from Noam)
end
t_read_database = cputime - t_read_database

num_species = size(Header_uni,1); % ~450,000 the entire database
if(plot_summary_statistics) % Just use some summary statistics (over the whole database
    plot_bacteria_summary_statistics(Header_uni, len_uni, Sequence_packed, bcs_figs_path)
end

if(num_species_in_database < num_species) % reduce database for faster reconstruction
    P = randperm(num_species);
    Header_uni = Header_uni(P(1:num_species_in_database));
    len_uni = len_uni(P(1:num_species_in_database));
    Sequence_packed = Sequence_packed(P(1:num_species_in_database));
    num_species = num_species_in_database;
end
clear posMix

[species_present I J] = intersect(Header_uni, Header_mix); posMix(J) = I; % Find where mixture baceria are in big database
posMix = posMix(posMix > 0); % Temp!!! (since we lost some species)
correctWeight = zeros(num_species,1);
correctWeight(posMix) = w(J)./sum(w(J)); % this is the solution: only 100 values > 0 (we may have lost some of them)



t_build_mixture_matrix = cputime;
[ReadBySpeciesMat values] = BuildMixingMatrixFromSequences(...
    readLength, Sequence_packed, len_uni, matlab_word_size); % find for last 1000
t_build_mixture_matrix = cputime - t_build_mixture_matrix

for i=1:length(algorithm_str_vec) % try all algorithms 
    run_alg = algorithm_str_vec{i}
    ReconstructMixtureFromReads(0, ReadBySpeciesMat, values, ...
        uniqueReads, uniqueReads_inds, 1); % This is the main 'master' function 
end



k = 1;
curr_kp = 1:num_species; % vector of bacterias are left
while (length(curr_kp)>1000) % what's going on here?
    currX = iterate(uniqueReads, uniqueReads_counts, ...
        Sequence_packed(curr_kp), Header_uni(curr_kp), readLength, ...
         algorithm_str, block_size, []);
    next_kp = curr_kp(find(currX>10^-3)); % find all variables larger than a pre-defined threshold
    store_kp{k} = curr_kp;
    store_X{k} = currX;
    k = k +1;
    curr_kp = next_kp; % hold indices of the remaining species (less and less)
    save(fullfile(bcs_output_path, 'iterations100outof450000_based10e6_similar'), ...
        'store_kp', 'store_X', 'posMix', 'correctWeight');  % save intermediate results
end
store_kp{k} = curr_kp;
store_X{k} = [];
save(fullfile(bcs_output_path, 'iterations100outof450000_based10e6_similar'), ...
    'store_kp', 'store_X', 'posMix', 'correctWeight');  % save final results

% solve again with only the set of remaining bacteria
tmpInd = store_kp{length(store_kp)};
Sequence1 = Sequence_uni(tmpInd);
Header1 = Header_uni(tmpInd);
[ReadBySpeciesMat values] = BuildMixingMatrixFromSequences(readLength,Sequence1); % find for last 1000
z = SolveMixingMatrixFromReads(0,ReadBySpeciesMat,values,uniqueReads,uniqueReads_inds,1);
reconstructed_weights = zeros(1,num_species); reconstructed_weights(tmpInd) = z;

CompareSolutionToTrueMixture(reconstructed_weights, correct_weights, comparison_output_file);  % Compare original and resulting results

