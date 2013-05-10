function pos=get_mean_posalign(read1, database1,database_size)
% database = Sequence_uni_amp_uni_char_array(:); so it will be a one col
% char array
% database_size = size(Sequence_uni_amp_uni_char_array); this is the size
% of the array we need it for the ind2sub function
% the output is the mean postion for those who has alignment, if there is
% no alignment for the (+) strand it search for the reverse complement and
% give the mean position but +1000, so the output pos can be uint16


k=strfind(database1',read1);
if ~isempty(k)
    [I,~]=ind2sub(database_size,k);
    pos=mean(I);
else
    k=strfind(database1',seqrcomplement(read1));
    if ~isempty(k)
        [I,~]=ind2sub(database_size,k);
        pos=1000+mean(I);
    else
        pos=0;
    end    
end



