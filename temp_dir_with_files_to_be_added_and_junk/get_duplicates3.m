% the same as get_duplicates2 but outputs the required for myUniqueUINT8. Does not output the index

function [out] = get_duplicates3(v)
%keyboard

origLenV = length(v);
[v sort_perm] = sort(v);
[vals I J] = unique(v);

if length(vals)<origLenV
  %keyboard
  num_unique = length(vals);
  num_dups = myHistForBAC(J,1:num_unique);
else
  num_dups = ones(1,origLenV);
end

out = [vals';num_dups];
