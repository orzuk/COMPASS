function [pos sd rev_pos rev_sd]=get_mean_posalign_generic(read,database1,database_size,bacseq_len,searchForwardAndReverseFlag)
% database1 = Sequence_uni_amp_uni_char_array(:); so it will be a one col
% char array so in order that work you need to make sure that
% Sequence_uni_amp_uni_char_array has the dimensions of seq_lengthXnumber_of
%_bacteria
% database_size = size(Sequence_uni_amp_uni_char_array); this is the size
% of the array we need it for the ind2sub function
% the output is the mean postion for those who has alignment, if there is
% no alignment for the (+) strand it search for the reverse complement and
% give the mean position but +1000, so the output pos can be uint16



if searchForwardAndReverseFlag==1 % DNA
  k=strfind(database1',read);
  if ~isempty(k)
    [I,~]=ind2sub(database_size,k);
    pos=mean(I);
    sd = std(I);
    rev_pos = mean(I);
    rev_sd = sd;
    
  else
    k=strfind(database1',seqrcomplement(read));
    %keyboard
    if ~isempty(k)
      [I,J]=ind2sub(database_size,k);
      pos=10000+mean(bacseq_len(J)-I);
      sd = 10000+std(bacseq_len(J)-I);
      
      rev_pos = 10000+mean(I);
      rev_sd = 10000+std(I);
      
    else
      rev_pos = 0;
      rev_sd = 0;
      pos=0;
      sd = 0;
    end    
  end
else % RNA like
  k=strfind(database1',read);
  if ~isempty(k)
    [I,~]=ind2sub(database_size,k);
    pos=mean(I);
    sd = std(I);
    rev_pos = mean(I);
    rev_sd = sd;
  else
    rev_pos = 0;
    rev_sd = 0;
    pos=0;
    sd = 0;
  end
end








