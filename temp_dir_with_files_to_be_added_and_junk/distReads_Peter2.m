function [u_right,b]=distReads_Peter2(POS,uniqueReads_length)

POS = round(POS);

a = find(POS);

b = unique(POS(a));
num_right = zeros(length(b),1);
for i=1:length(b)
  c = find(POS(a)==b(i));
  num_right(i) = sum(uniqueReads_length(a(c)));
end



