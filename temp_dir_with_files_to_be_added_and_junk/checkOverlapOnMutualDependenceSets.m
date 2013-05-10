function new_s=checkOverlapOnMutualDependenceSets(region_read,region,redM,s1,s2)

%keyboard

[x,y_s1] = find(region(s1,:));y_s1 = unique(y_s1);
[x,y_s2] = find(region(s2,:));y_s2 = unique(y_s2);

region_inter = intersect(y_s1,y_s2);

r = [];
for i=1:length(region_inter)
  r = [r;region_read{region_inter(i)}];
end

new_s = findLinearDependenceSets(redM(r,[s1,s2]));

orig = [s1,s2];
for i=1:length(new_s)
  new_s{i} = orig(new_s{i});
end

