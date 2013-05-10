% Simulate reads: (output is unique reads and their counts) 
function [uniqueReads,uniqueReads_length]=createReads_package(auxData)

[simulated_reads] = createReadsForSpecificMixture_package(auxData);

% create a list of unique reads (why not just perform unique?)
[uniqueReads,uniqueReads_inds] = extract_sub_kmers(simulated_reads, auxData.readLength*ones(size(simulated_reads,1),1),auxData.readLength, 1,0);
clear simulated_reads

% find the number of appearances of each read
[junk_vals junk_inds uniqueReads_length] = get_duplicates(uniqueReads_inds(:,1));% uniqueReads_length is number of appearances of each read


