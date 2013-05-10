function basic=listAmbiguousSeq(currSeq)
%keyboard
amb = seq2regexp(currSeq,'ambiguous','false');

a = findstr(amb,'[');
b = findstr(amb,']');

a(end+1) = length(amb)+1;
b(end+1) = length(amb)+1;

%group = cell(1,length(a));
basic = currSeq(1:a(1)-1);
for i=1:length(a)-1
  group = amb(a(i)+1:b(i)-1);
  
  tmpBasic = repmat(basic,length(group),1);
  
  add = repmat(group,size(tmpBasic,1)/length(group),1);
  
  tmpBasic(:,end+1) = add(:)';
  
  next = amb(b(i)+1:a(i+1)-1);
  tmpNext = repmat(next,size(tmpBasic,1),1);
  basic = [tmpBasic,tmpNext];
  %pause
end



%size(basic)
%size(currSeq)
