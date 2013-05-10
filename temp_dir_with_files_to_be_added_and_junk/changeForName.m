function out=changeForName(in)

out = [];
for i=1:length(in)
  out = [out,'_',num2str(in(i))];
end
