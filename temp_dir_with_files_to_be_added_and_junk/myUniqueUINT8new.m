function [vals,inds,positions,leng]=myUniqueUINT8new(seqO,posO)
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


b = find(c~=0);
clear c  
%keyboard
b = [0;b];
if b(end)~=length(index)
  b = [b;length(index)];
end

inds = cell(length(b)-1,1);
vals = zeros(length(b)-1,size(seqO,2),'uint8');


%keyboard
if exist('posO')
  
  positions = inds;
  
  oneInB = find(diff(b)==1)';
  other = 1:length(b)-1;
  other(oneInB) = [];
  
  leng = diff(b);
  inds = mat2cell(index,leng,1);
  vals = seqO(index(b(2:end)),:);

  lenOneInB = length(oneInB);
  positions(oneInB) = mat2cell([posO(index(b(oneInB+1)))';ones(1,lenOneInB)],2,ones(1,lenOneInB));
  %keyboard
  
  %z = [b(other)+1;b(other+1)];
  %tryPos = mat2cell(posO(index(z(:,1):z(:,2))),leng((other)),1);
  %positions(other) = cellfun(@get_duplicates3,tryPos,'UniformOutput',0);
  
  tryPos = mat2cell(posO(index),leng,1);
  positions(other) = cellfun(@get_duplicates3,tryPos(other),'UniformOutput',0);
 

  
  
else % single output
  positions = [];
  leng = [];
  for i=1:length(b)-1
    inds{i} = index(b(i)+1:b(i+1));
    vals(i,:) = seqO(index(b(i)+1),:);
  end

end


