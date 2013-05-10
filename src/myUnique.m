% Perform unique but keep indices
%
% Input:
% seqO - vector of sequences 
% posO - vector of positions 
%
% Output:
% vals - unique values  
% inds - indices of unique values 
% positions - where ..
% leng - ?? 
%
function [vals,inds,positions,leng]=myUnique(seqO,posO)

[junk,index] = sortrows(seqO);
a = abs(diff(abs(junk)));
clear junk
c = sum(a,2);
clear a
b = find(c~=0);
clear c
%keyboard
b = [0;b];
if b(end)~=length(index)
    b = [b;length(index)];
end
inds = cell(length(b)-1,1);
vals = char(55*ones(length(b)-1,size(seqO,2))); % this is something good only for chars 

if exist('posO', 'var') % also consider positions ?? 
    positions = inds;
    leng = zeros(length(b)-1,1);
    for i=1:length(b)-1
        inds{i} = index(b(i)+1:b(i+1));
        [vals1 inds1 num_dups] = get_duplicates(posO(inds{i}));
        positions{i} = [vals1';num_dups];
        leng(i) = length(vals1);
        vals(i,:) = seqO(index(b(i)+1),:);
    end
else
    positions = [];
    leng = [];
    for i=1:length(b)-1
        inds{i} = index(b(i)+1:b(i+1));
        vals(i,:) = seqO(index(b(i)+1),:);
    end    
end

