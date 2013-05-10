function [u_right,n_right]=distReads_Peter(right,right_res)

a = find(right_res.POS);
length(a)


b = unique(right_res.POS(a));
num_right = zeros(length(b),1);
for i=1:length(b)
  c = find(right_res.POS(a)==b(i));
  num_right(i) = sum(right.uniqueReads_length(a(c)));
end

%keyboard
z_right = round(b);
u_right = unique(z_right);
for i=1:length(u_right)
  n_right(i) = sum(num_right(find(z_right==u_right(i))));
end
