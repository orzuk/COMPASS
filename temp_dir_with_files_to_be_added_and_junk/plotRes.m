function plotRes(name,resName)

origName = name;
load(name)

a = find(name=='/');
name(1:a(end)) = [];
name(find(name=='_')) = '-';
load(resName)
found
%keyboard
for i=1:length(found)
  if ~isempty(found{i})
    z = find(found{i}<0);
    figure(i)
    plot(abs(found{i}),correctWeight,'.');
    hold on
    plot(abs(found{i}(z)),correctWeight(z),'ro');
    title([name,num2str(i)])
  end
end
%keyboard
%find(found{4}<1e-4 & correctWeight'>0.005)
