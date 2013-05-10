% Debug matrix building: look at Noam's and my version and see which one is better
% function debug_build_matrix()
clear

AssignGeneralConstants;

%userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';

switch machine
    case PC
        userDir = '../../compressed_sensing/metagenomics/next_gen';
    case UNIX
        userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
end
addpath([userDir,'/CS/mFilesBAC'])


basicSeqKey = [userDir,'/CS/BAC/dat450000/key450000'];
basicSeqNameDir = [userDir,'/CS/BAC/dat450000/'];
readLength = 50;

%b = 1:5; % 1674; % 1:5;
b = 1674; % 1:5;
load10=0;
if(load10)
    load('seq1.mat'); % just first 10 
    Sequence_packed = pack_seqs(Sequence1, matlab_word_size);
else
%    load('first2000.mat'); % first 2000
    load('first10000.mat'); % first 10000
    Sequence1 = Sequence_uni; 
    Sequence_packed = pack_seqs(Sequence1, matlab_word_size);
    
%     seqs_len = (64/2)*length_cell(Sequence_uni);
%     Sequence_packed = Sequence_uni;
%     for i=1:length(Sequence_packed)
%         Sequence1{i} = unpack_seqs(Sequence_packed{i}, seqs_len(i), 64);
%     end
end
len_uni = length_cell(Sequence1);
tic
[OrMat values] = BuildMixingMatrixFromSequences(readLength, Sequence_packed(b), len_uni(b), matlab_word_size); % find for last 1000
OrKmers = int2nt(unpack_seqs(values, readLength, matlab_word_size));
toc



tic
[NoamMat NoamKmers] = prepareGroupOf1000DistributedSequenceFiles2(readLength,b,basicSeqNameDir,basicSeqKey);
NoamKmers = char(NoamKmers);
for i=length(b)
    NoamMat(:,i) = round(NoamMat(:,i) .* (len_uni(b(i))-readLength+1)); % get back integers 
end
toc



%load(fullfile(user_dir, '~/CS/BAC/s16_data_uni_packed')
% if(~exist('Sequence_packed', 'var'))
%     load([userDir,'/data/s16_data_uni_packed']);
% end
% Sequence_packed1 = Sequence_packed(b);


for i=b
    s{i} = unpack_seqs(Sequence_packed{i},len_uni(i), matlab_word_size);
end

a = zeros(size(NoamMat));
a(find(NoamMat)) =1;
z = sum(a,2);

a = zeros(size(OrMat));
a(find(OrMat)) =1;
w = sum(a,2);

q = unpack_seqs(values(50,:),readLength, matlab_word_size);
q = int2nt(q);
[junk,i1,i2]=intersect(NoamKmers,q,'rows');
size(NoamMat)
size(OrMat)


% Align kmers
[NoamKmersSorted NoamPerm] = sortrows(NoamKmers);
[OrKmersSorted OrPerm] = sortrows(OrKmers);
%OrKmersSorted(4,5) = 'F';

only_or_kmers = setdiff(OrKmers, NoamKmers, 'rows')
only_noam_kmers = setdiff(NoamKmers, OrKmers, 'rows')
different_kmers = setxor(OrKmers, NoamKmers, 'rows');

OrCounts = sum(OrMat(OrPerm,:),2);
NoamCounts = sum(NoamMat(NoamPerm,:),2);
different_kmers_sorted = find(NoamKmersSorted ~= OrKmersSorted, 1)
if(machine == PC)
    figure; imagesc(OrMat); colorbar;
    figure; imagesc(OrMat(OrPerm,:) - NoamMat(NoamPerm,:)); colorbar; 
    title('Difference between two methods');     % should be all zero !!!
end
max_error = max(max(abs(OrMat(OrPerm,:) - NoamMat(NoamPerm,:))))


total_or_counts = sum(OrCounts)
total_noam_counts = sum(NoamCounts)

mismatch_kmers_inds = find(OrCounts ~= NoamCounts);
