% FindSeqInDB
% find a close sequence in the database to a given sequence.
% requires a prepared local blast database named 454idx.fa
% input
% seq - the sequence to look for
% output
% idx - the index (name in the database) of the close sequence
% or 0 if no similar sequence (>=98% identity) is found

function [idx]=myFindSeqInDB_using_database_with_ambiguous(seq,basicDir,i,k)



filename=[basicDir,'/tmp-','_',num2str(i),'_',num2str(k),'_',num2str(randi(1E8)) '.fa'];
fastawrite(filename,'head',seq);

results = blastlocal('InputQuery',filename, 'Program', 'blastn','Database','/homes/csfaculty/shental/CS/BAC/12Samples/454/idx/All454Idx.fa');

delete(filename);
%keyboard
ba=results.Hits(2).HSPs.Alignment;
periden=sum(ba(2,:)=='|')/length(seq);
%keyboard
if (periden>=0.98)
    idx=str2num(results.Hits(1).Name)
else
    idx=0;
end

