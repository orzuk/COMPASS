function F=findSpecificBacteria(a,Header_no_nonACGT)

%keyboard
F = zeros(length(a),1);
for i=1:length(a)
  i
  flag = 0;
  k = 1;
  while flag==0 && k<=length(Header_no_nonACGT)
    for j=1:length(Header_no_nonACGT{k})
      if ~isempty(findstr(Header_no_nonACGT{k}{j},a{i}))
        F(i) = [k];
        flag = 1;
      end
    end
    k = k+1;
  end
end
