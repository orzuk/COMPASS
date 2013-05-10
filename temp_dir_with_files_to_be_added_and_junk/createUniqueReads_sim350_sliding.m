function createUniqueReads_sim350_sliding(dirName,fileName,readLength)

load([dirName,'/',fileName]);

%keyboard

if ~isempty(find(size(reads_uni)==1))
  uniqueReads_length = reads_uni;
  uniqueReads = reads_uni_freq;
else
  uniqueReads_length = reads_uni_freq;
  uniqueReads = reads_uni;
end

uniqueReads = pack_seqs(uniqueReads,64);
whos
%keyboard
save([dirName,'/reads/',fileName,'_sliding',num2str(readLength)],'uniqueReads','uniqueReads_length');
pause(2)

%keyboard