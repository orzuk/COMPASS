function [vals,inds,positions,leng]=myUniqueUINT8(seqO,posO)
% vals = the unique rows in Esq
% currents = the indices of the unique rows in seqO
% positions = the BAC which holds the sequence and the number of times it does.
% leng = the total number of BAC which hold the sequence

disp('assumes that the first input to be sorted is uint8. myUniqueUINT8.m')
%keyboard
%[junk1,index1] = sortrows(seqO);
[index] = sortrowsc(seqO,1:size(seqO,2));

disp('myUniqueUINT8 - maybe start with a smaller pt is seqO is too large. myUniqueUINT8.m')
okFlag = 0;
pt = size(seqO,1)-1;
while okFlag==0
  try
    part = 1:pt:size(seqO,1);
    part(end) = size(seqO,1)+1;
    
    c = zeros(size(seqO,1)-1,1);
    for i=1:length(part)-1
      
      tmpInd = [part(i):part(i+1)-1]';
      a = abs(diff(abs(seqO(index(tmpInd),:))));
      c(tmpInd(1:end-1)) = sum(a,2);
      if i>1
        a = abs(diff(abs(seqO(index(tmpInd-1:tmpInd),:))));
        c(tmpInd(1)-1) = sum(a,2);
      end
    end
    okFlag = 1;
  catch me
    pt = round(pt/2)
    okFlag = 0;
  end
end

clear a 


%a = abs(diff(abs(junk)));
%clear junk
%c = sum(a,2);
%clear a
b = find(c~=0);
clear c  
%keyboard
b = [0;b];
if b(end)~=length(index)
  b = [b;length(index)];
end

inds = cell(length(b)-1,1);
%vals = char(55*ones(length(b)-1,size(seqO,2)));
vals = zeros(length(b)-1,size(seqO,2),'uint8');


%keyboard
if exist('posO')
  positions = inds;
  leng = zeros(length(b)-1,1);
  
  % single
  oneInB = find(diff(b)==1)';
  leng(oneInB) = 1;
  inds(oneInB) = num2cell(index(b(oneInB+1)));
  
  lenOneInB = length(oneInB);
  positions(oneInB) = mat2cell([posO(index(b(oneInB+1)))';ones(1,lenOneInB)],2,ones(1,lenOneInB));
  vals(oneInB,:) = seqO(index(b(oneInB)+1),:);
  
  
    
  % several
  other = 1:length(b)-1;
  other(oneInB) = [];
  for i=other
    inds{i} = index(b(i)+1:b(i+1));
    [vals1 inds1 num_dups] = get_duplicates2(posO(inds{i}));    
    positions{i} = [vals1';num_dups];
    leng(i) = length(vals1);
    vals(i,:) = seqO(index(b(i)+1),:);
  end
  
else % single output
  positions = [];
  leng = [];
  for i=1:length(b)-1
    inds{i} = index(b(i)+1:b(i+1));
    vals(i,:) = seqO(index(b(i)+1),:);
  end

end


