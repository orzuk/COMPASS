function s=findLinearDependenceSets(redM)
%keyboard
B=full(redM'*redM);
D=rref_mod2(B);

k = 0;
p = zeros(length(find(D)),2);
for i=1:size(D,1)
  a = find(D(i,:));
  if ~isempty(a)
    if length(a)==1
      add_p = [a,a];
    else
      add_p = nchoosek(a,2);
    end
    l = size(add_p,1);
    p(k+1:k+l,:) = add_p;
    k = k+l;
  end
end
p(k+1:end,:) = [];

a = unionfind(p,size(D,1));
%keyboard
u = unique(a);
s = cell(1,length(u));
for i=1:length(u)
  s{i} = find(a==u(i));
end

%keyboard

rev_s = [];
for i=1:length(u)
  %i
  out = findPositiveLinDep(s{i},redM);
  rev_s = [rev_s,out];
end

s = rev_s;