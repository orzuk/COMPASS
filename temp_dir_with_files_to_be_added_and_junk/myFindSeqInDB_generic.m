% FindSeqInDB
% find a close sequence in the database to a given sequence.
% requires a prepared local blast database named 454idx.fa
% input
% seq - the sequence to look for
% output
% idx - the index (name in the database) of the close sequence
% or 0 if no similar sequence (>=98% identity) is found

function [idx]=myFindSeqInDB_generic(basicName,basicDir,databaseName,blastPath,indexFirst,indexLast,saveFileName,dataDir,dataName)

%rmpath(genpath('/home/csfaculty/shental/Weizmann/ph2users/fenoam/'));rmpath('/home/csfaculty/shental/matlab');cd ~/CS/mFilesBAC;addpath('~/CS/mFilesBAC');  mcc -v -m myFindSeqInDB_generic.m  -d /homes/csfaculty/shental/CS/BAC/cvx 

% databaseName = [userDir,'~/CS/BAC/12Samples/454/idx/old_fasta_files_without_ambiguous/454idx.fa'];

if isdeployed%ischar(indexFirst)
  indexFirst = str2num(indexFirst); 
  indexLast = str2num(indexLast);
end

% load reads
% load([dataDir,'/',dataName],'reads_uni');

% change 20.6
load([dataDir,'/',dataName]);


if ~exist('lenSeq')
  lenSeq = size(reads_uni,2)*ones(1,size(reads_uni,1));
end
%keyboard
% go over the list - with path
PWD = pwd;
cd(blastPath)
for i=indexFirst:indexLast
  idx(i) = findIn(basicDir,basicName,databaseName,blastPath,reads_uni(i,1:lenSeq(i)),i);
end
cd(pwd);

save(saveFileName,'idx','indexFirst','indexLast')

function idx=findIn(basicDir,basicName,databaseName,blastPath,seq,i)
%keyboard
filename=[basicDir,'/tmp/tmp_',basicName,'_',num2str(i),'_',num2str(randi(1E8)) '.fa'];
fastawrite(filename,'head',seq);


results = blastlocal('InputQuery',filename, 'Program', 'blastn','Database',databaseName,'blastpath',[blastPath,'blastall']);


delete(filename);
%keyboard
ba=results.Hits(2).HSPs.Alignment;
periden=sum(ba(2,:)=='|')/length(seq);

if (periden>=0.98)
    idx=str2num(results.Hits(1).Name)
else
    idx=0;
end


