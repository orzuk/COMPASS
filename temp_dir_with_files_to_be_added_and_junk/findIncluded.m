function OneRegionIncludedInOther=findIncluded(new_Ap,regionForRead,numRegions)

% swollow those that are amplified only in 1 region
keyboard
singleRegBact = find(sum(new_Ap,1)==2*numRegions); 
for i=singleRegBact
  p = find(new_Ap(:,i));
  if length(p)~=2*numRegions
    'problem'
    pause
  end
  if length(unique(regionForRead(p)))~=numRegions
    'problem'
    pause
  end
end

% 

clear other independent
independent = [];
OneRegionIncludedInOther = cell(1);
k = 1;
%keyboard
for i=singleRegBact
  x = find(new_Ap(:,i));
  
  % is it included in any other? namley both forward and reverse are part of another
  other1 = find(new_Ap(x(1),1:end-6));
  for j=2:length(x)
    other1 = intersect(other1,find(new_Ap(x(j),1:end-6)));
  end
  
  tmpOther = setdiff(other1,i);
  if isempty(tmpOther) % not included in any other - check if stands by itself
    val = full(abs(new_Ap(x,size(new_Ap,2)-6+regionForRead(x(1)))));
    %val
    independent = [independent,i];
  else
    %1
    %pause
    OneRegionIncludedInOther{k} = [i,tmpOther];
    k = k+1;
  end
%  pause
end

