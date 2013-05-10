function name=addNamesForMichael_in(nm,Header_750)

name = cell(1,length(nm));
curr_nm = nm;
for i=1:length(Header_750)
  if mod(i,100)==0
    i
  end
  if length(curr_nm)>0
    for j=1:length(Header_750{i})
      currNumber=findNumberFromHeader(Header_750{i});
      [junk,i1,i2] = intersect(nm,currNumber);
      if ~isempty(junk)
        name{i1} = Header_750{i};
        curr_nm(find(curr_nm==currNumber)) = [];
      end
    end
  end
end



