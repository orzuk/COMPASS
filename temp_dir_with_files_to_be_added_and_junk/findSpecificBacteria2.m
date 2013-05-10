function found=findSpecificBacteria2(Header_uni,nm)

%keyboard
l = 1;
for i=1:length(Header_uni)
  for j=1:length(Header_uni{i})
    if iscell(Header_uni{i}{j})
      for k=1:length(Header_uni{i}{j})
        f = findstr(Header_uni{i}{j}{k},nm);
        if ~isempty(f)
          Header_uni{i}{j}{k}
          found(l) = [i];return
          l = l+1;
          %pause
        end
      end
    else
      f = findstr(Header_uni{i}{j},nm);
      if ~isempty(f)
        Header_uni{i}{j}
        found(l) = [i];return
        l = l+1;
        %pause
      end
    end
    
    
  end
end
