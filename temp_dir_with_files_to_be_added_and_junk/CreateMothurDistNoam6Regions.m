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

function [distmat,allSeqs,origFreqVec,recFreqVecL1,recFreqVecL2]=CreateMothurDistNoam6Regions(fileName,found,correctWeight,auxData)

%keyboard


allSeqs = unique([find(correctWeight),find(found)]);



% Create the joined vector and the appropriate frequency vectors

% load the sequences
disp(['reading ' num2str(length(allSeqs)) ' sequences']);
basicSeqNameDir = auxData.basicSeqNameDir;
basicSeqKey= auxData.basicSeqKey;

a = findstr(fileName,'/');
tmpfilename = ['seqs_',fileName(a(end)+1:end)];

load(basicSeqKey,'len_uni')

% remove the fastawrite append warnings
warnState = warning; %Save the current warning state
warning('off','Bioinfo:fastawrite:AppendToFile');

%load the sequences and save to fasta file
delete(tmpfilename);
tmpInd = allSeqs;
[HeaderAll,SequenceAll] = loadSeq(tmpInd,basicSeqNameDir,basicSeqKey);
for a=1:length(tmpInd)
    allCSeq=int2nt(unpack_seqs(SequenceAll{a},len_uni(tmpInd(a)) , 64));
    fastawrite([tmpfilename,'.fa'],['Seq-' num2str(a)],allCSeq);
end

% and restore warning state
warning(warnState); %Reset warning state to previous settings

% run the mothur for multiple sequence alignment and distance matrix
% calculation
disp('starting mothur');
%keyboard

pwd1 = pwd;
cd(auxData.currDir);

w = ['!',auxData.userDir,'/CS/BAC/mothur/mothur "#align.seqs(candidate=',...
     tmpfilename,'.fa'...
     ', template=',auxData.userDir,'/CS/BAC/core_set_aligned.fasta.imputed); dist.seqs(fasta=',...
     tmpfilename,'.align, output=square)"'];
eval(w)

disp('Finished calculating distances');

% load the mothur output matrix and skip the first line
distmat=dlmread([auxData.currDir,'/',tmpfilename,'.square.dist'],'\t',0,1);
distmat=distmat(2:end,1:end-1);

% save the output
%save([fileName,'_mothurRes'],'distmat','allSeqs','origFreqVec','recFreqVec');
save([fileName,'_mothurResL1L2'],'distmat','allSeqs');



unix(['cd ',auxData.currDir,';rm -f ',tmpfilename,'.fa ',tmpfilename,'.align.report ',tmpfilename,'.align ',tmpfilename,'.square.dist mothur.*.logfile'])
cd(pwd1);
