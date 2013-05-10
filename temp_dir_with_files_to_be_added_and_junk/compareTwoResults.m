function [distmat,allSeqs,origFreqVec,recFreqVec]=compareTwoResults(fileName1,fileName2,minFreqToConsider,auxData)


threshold=0;


r1 = load([fileName1]);
r2 = load([fileName2]);

disp('compares the last!!! maybe change. CreateMothurDistTwoSolutions.m')
orig = find(abs(r1.found{end})>minFreqToConsider);
rec = find(abs(r2.found{end})>minFreqToConsider);
%keyboard


% load the sequences
basicSeqNameDir = auxData.basicSeqNameDir;
basicSeqKey= auxData.basicSeqKey;
load(basicSeqKey,'len_uni')


[HeaderAll,SequenceOrig] = loadSeq(orig,basicSeqNameDir,basicSeqKey);
[HeaderAll,SequenceRec] = loadSeq(rec,basicSeqNameDir,basicSeqKey);
%keyboard
ham = NaN*ones(length(orig),length(rec));
for i=1:length(orig)
  seq_orig = int2nt(unpack_seqs(SequenceOrig{i},len_uni(orig(i)),64));
  for j=1:length(rec)
    seq_rec = int2nt(unpack_seqs(SequenceRec{j},len_uni(rec(j)),64));
    [cscore,algn]=swalign(seq_orig,seq_rec,'Alphabet','NT');
    ham(i,j) = length(find(algn=='|'));
  end
end

imagesc(ham)
xlabel(fileName1);
ylabel(fileName2);
