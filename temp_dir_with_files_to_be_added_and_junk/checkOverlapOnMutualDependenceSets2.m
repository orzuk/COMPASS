function new_s=checkOverlapOnMutualDependenceSets2(region_read,region,redM,varargin)
% the same as checkOverlapOnMutualDependenceSets but checks sets of groups
%keyboard

region_inter = 1:6;
bact = [];
for i=1:length(varargin{1})
  [x,y] = find(region(varargin{1}{i},:));y = unique(y);
  region_inter = intersect(region_inter,y);
  bact = [bact,varargin{1}{i}];
end
bact = unique(bact);

r = [];
for i=1:length(region_inter)
  r = [r;region_read{region_inter(i)}];
end

new_s = findLinearDependenceSets(redM(r,bact));


for i=1:length(new_s)
  new_s{i} = bact(new_s{i});
end

