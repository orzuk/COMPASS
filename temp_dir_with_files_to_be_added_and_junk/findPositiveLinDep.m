function new_s=findPositiveLinDep(s,redM)
%keyboard
z = redM(:,s); D = rref_mod2(z'*z);
[x,y] = find(z);

zz = full(z(x,:));

res = [1:size(zz,2)]';
for k=1:size(zz,2)
  rest = 1:size(zz,2);
  rest(k) = [];
  curr_s = linsolve(zz(:,rest),zz(:,k));
  if isempty(find(curr_s<-10^-3));
    out = [k,rest(find(curr_s'>10^-3))];
    a = min(res(out));
    res([out find(res==a)']) = a;
  end
end

[junk,new_s,junk2] = get_duplicates(res);
new_s = cellfun(@(x) s(x),new_s,'uniformoutput',false);
new_s = new_s';
