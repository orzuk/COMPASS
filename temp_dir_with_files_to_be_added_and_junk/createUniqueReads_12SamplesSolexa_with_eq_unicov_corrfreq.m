
function createUniqueReads_12SamplesSolexa_with_eq_unicov_corrfreq(dirName,readLength,samplesName,typeName)

load(['~/CS/BAC/12Samples/Solexa/data/',dirName,'/illumina_reads',num2str(readLength),'_',samplesName,typeName],'uni_read_seq','freq_read','freq_corr2pos');

%keyboard

%size(uni_read_seq)
%[u i1 i2]= unique(uni_read_seq,'rows');
%size(u)
%disp('%%%')
%return


% see if there are non ACGT
if 1==2
  n = zeros(size(uni_read_seq,1),1);
  for i=1:size(uni_read_seq,1)
    if ~isempty(find(uni_read_seq(i,:)~='A' & uni_read_seq(i,:)~='C' & uni_read_seq(i,:)~='G' & uni_read_seq(i,:)~='T'))
      n(i) = 1;
    end
  end
  length(find(n))
  return
end

uniqueReads = pack_seqs(uni_read_seq,64);

% with equalization
uniqueReads_length = freq_corr2pos';
save(['~/CS/BAC/12Samples/Solexa/data/',dirName,'/readsSolexa_with_equalization_',dirName,'_',num2str(readLength),'_',samplesName],'uniqueReads','uniqueReads_length');

%keyboard
