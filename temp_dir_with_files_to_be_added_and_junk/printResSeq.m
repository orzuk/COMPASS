function printResSeq(b,results,Header_uni,fileName)
%keyboard
fid = fopen(fileName,'w');
for i=1:length(b)
  
  i
  
  nm = [];
  for j=1:length(Header_uni{b(i)})
    if iscell(Header_uni{b(i)}{j})
      for k=1:length(Header_uni{b(i)}{j})
        nm = [nm,';',Header_uni{b(i)}{j}{k}];
        %Header_uni{b(i)}{j}{k}
      end
    else
      nm = [nm,';',Header_uni{b(i)}{j}];
    end
  end
  %nm
  fprintf(fid,'%s \t %1.2f\n',nm,results(i));
end
fclose(fid);