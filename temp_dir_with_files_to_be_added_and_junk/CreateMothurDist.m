% calculate the pairwise distance matrix using the mothus program
% input:
% dirnane,filename - where the algorithm output is stored
% outfilename - filename where to store the distance matrix
% israndom - 0 do normal analysis, 1 randomize the recreated vector
%
% output:
% distmat - the pairwise mothur distance matrix
% allSeqs - vector containing the indices of the sequences of the distance
% matrix (union of the original and recreated vectors)
% origFreqVec - vector containing the frequencies of the original sequences
% in the distance matrix
% recFreqVec - vector containing the frequencies of the recreated sequences
% in the distance matrix

function [distmat,allSeqs,origFreqVec,recFreqVec]=CreateMothurDist(dirName,fileName,outfilename,israndom)

threshold=0;

r = load([dirName fileName]);

orig=r.resCell{1};
rec=r.resCell{2};

origSeqs=orig(:,1);
origFreq=abs(orig(:,2));
recSeqs=rec(:,1);

    %RANDOM: REMOVE IF WANT NON RANDOM!!!!!!!!!!!!!
if (israndom==1)
    disp('Random iteration!');
    recSeqs=randi(410849,length(recSeqs),1);
end

recFreq=abs(rec(:,2));
% remove low abundance sequences (see threshold value)
recSeqs=recSeqs(recFreq>threshold);
recFreq=recFreq(recFreq>threshold);

% Create the joined vector and the appropriate frequency vectors
allSeqs=union(origSeqs,recSeqs);
origFreqVec=zeros(size(allSeqs));
recFreqVec=zeros(size(allSeqs));

for a=1:length(origSeqs)
    origFreqVec(find(allSeqs==origSeqs(a)))=origFreq(a);
end
for a=1:length(recSeqs)
    recFreqVec(find(allSeqs==recSeqs(a)))=recFreq(a);
end

% load the sequences
disp(['reading ' num2str(length(allSeqs)) ' sequences']);
basicSeqNameDir = ['..\fig2\packed64\'];
basicSeqKey= ['..\fig2\keyNoNonACGT'];

tmpfilename='seqs.fa';

load(basicSeqKey,'len_uni')

% remove the fastawrite append warnings
warnState = warning; %Save the current warning state
warning('off','Bioinfo:fastawrite:AppendToFile');

%load the sequences and save to fasta file
delete(tmpfilename);
tmpInd = allSeqs;
[HeaderAll,SequenceAll] = loadSeqNames(tmpInd,basicSeqNameDir,basicSeqKey);
for a=1:length(tmpInd)
    allCSeq=int2nt(unpack_seqs(SequenceAll{a},len_uni(tmpInd(a)) , 64));
    fastawrite(tmpfilename,['Seq-' num2str(a)],allCSeq);
end

% and restore warning state
warning(warnState); %Reset warning state to previous settings

% run the mothur for multiple sequence alignment and distance matrix
% calculation
disp('starting mothur');
!mothur\mothur "#align.seqs(candidate=seqs.fa, template=core_set_aligned.fasta.imputed); dist.seqs(fasta=seqs.align, output=square)"
disp('Finished calculating distances');

% load the mothur output matrix and skip the first line
distmat=dlmread('seqs.square.dist','\t',0,1);
distmat=distmat(2:end,1:end-1);

% save the output
save(outfilename,'distmat','allSeqs','origFreqVec','recFreqVec');
