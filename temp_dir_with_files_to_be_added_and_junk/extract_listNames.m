function list=extract_listNames(list50)

a = findstr(list50,'/seq/');
a = [a,length(list50)+1];
k = 1;
list = [];
for i=1:length(a)-1
  list{i} = [list50(a(i):a(i+1)-2)];
end
