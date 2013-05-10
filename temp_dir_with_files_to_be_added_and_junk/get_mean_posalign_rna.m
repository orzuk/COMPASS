function [pos sd]=get_mean_posalign_rna(read, database1,database_size,bacseq_len)
% database1 = Sequence_uni_amp_uni_char_array(:); so it will be a one col
% char array so in order that work you need to make sure that
% Sequence_uni_amp_uni_char_array has the dimensions of seq_lengthXnumber_of
%_bacteria
% database_size = size(Sequence_uni_amp_uni_char_array); this is the size
% of the array we need it for the ind2sub function
% the output is the mean postion for those who has alignment, if there is
% no alignment for the (+) strand it search for the reverse complement and
% give the mean position but +1000, so the output pos can be uint16

pos=0;
k=strfind(database1',read);
keyboard
if ~isempty(k)
    [I,~]=ind2sub(database_size,k);
    pos=mean(I);
    sd = std(I);
end


for i=1:size(Sequence_uni_chararray,2)
  if ~isempty(findstr(Sequence_uni_chararray(:,i)',read))
    i
  end
end
