function listThem(fid,sleepTime,list50)

a = findstr(list50,'/seq/');
a = [a,length(list50)+1];
k = 1;

for i=1:length(a)-1
  w = [list50(a(i):a(i+1)-2),' ;sleep 20;'];
  fprintf(fid,'%s\n',w);
  if mod(k,5)==0
    fprintf(fid,'sleep %s;\n',sleepTime);
  end
  k = k+1;
end
