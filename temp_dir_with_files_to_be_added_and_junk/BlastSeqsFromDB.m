% do a blast for sequences in the database and then find the index via the
% header
% input:
% seqs : a cell array of the sequences to search (text strings)
% dbheaders : the cell list of the database headers
% output:
% seqinds : indices of the sequences. 0 if not found
% matchscore : how good was the match from the db to the sequence
function [seqinds,matchscore]=BlastSeqsFromDB(seqs,dbheaders)
%% prepare the fasta file for the sequences we're searching for
FILENAME='tmpseqs.fasta';
LOCALDBNAME='seqs.fa'; % the blast database
keyboard
seqinds=zeros(length(seqs),1);
matchscore=zeros(length(seqs),1);

delete(FILENAME);
for a=1:length(seqs)
    cseq=seqs{a};
    chead=num2str(a);
    fastawrite(FILENAME,chead,cseq);
end

%% Do local blast to search for them
data=blastlocal('InputQuery', FILENAME,'Program','blastn','Database',LOCALDBNAME);
for a=1:length(data)
    cseqname=data(a).Hits(1).Name;
    pos=FindNameInDB(cseqname,dbheaders);
%    disp([cseqname num2str(pos)]);
    seqinds(a)=pos;
    matchscore(a)=data(a).Hits(1).HSPs(1).Identities.Percent;
end
