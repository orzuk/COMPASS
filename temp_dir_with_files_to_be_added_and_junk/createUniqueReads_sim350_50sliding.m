function createUniqueReads_sim350_50sliding(dirName)


load(['~/CS/BAC/sliding50_454/',dirName]);


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
save(['~/CS/BAC/sliding50_454/reads/',dirName,'_sliding50'],'uniqueReads','uniqueReads_length');
pause(2)

%keyboard