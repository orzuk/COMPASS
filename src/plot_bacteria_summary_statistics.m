% Just plot some simple statistics for bacterial database
%
% Input:
% species_names - names of bacterial species
% species_lengths - lengths of bacterial species
% packed_seqs - sequences themselfes
% bcs_output_path - where to save figures
%
function plot_bacteria_summary_statistics(species_names, species_lengths, ...
    packed_seqs, bcs_figs_path)

AssignGeneralConstants;
num_bins = 500;
[length_hist length_bin_locs] = hist(species_lengths, num_bins);
full_figure; bar(length_bin_locs, length_hist); % plot lengths distribution
mean_length = mean(species_lengths); median_length = median(species_lengths); std_lengths = std(species_lengths);

title(['Distribution of species lengths: mean=' ...
    num2str(mean_length, 4) ', median=' num2str(median_length, 4) ', std=' num2str(std_lengths,3)]);
xlabel('length (bp)'); ylabel('num. species');
my_saveas(gcf, fullfile(bcs_figs_path, 'bacterial_length_distribution'), format_fig_vec);




