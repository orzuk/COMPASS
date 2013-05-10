function plotGroup(group,len,name)

a = find(group);
if ~isempty(find(a==1))
  a(1) = [];
end

plot(len(a),group(a),'.');
xlabel('cluster size')
ylabel('number of relevant bacteria')
title([name,' #in not amp group: ',num2str(group(1))]);
