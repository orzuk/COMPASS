function printChromatogramResults(toPrint)
mx = cell(3,1);
for i=1:length(toPrint)
  for j=1:length(toPrint{i})
    mx{j} = max([mx{j},length(toPrint{i}{j})]);
  end
end

mx_t{1} = 10;
mx_t{2} = mx{1};
mx_t{3} = mx{2};

mx = mx_t;
s = length(toPrint);
for i=1:length(toPrint)
  for j=1:length(toPrint{i})
    a = toPrint{i}{j};
    switch a(1:2)
     case 'g_'
      a(1:2) = [];
      keyboard
      
      text(a,'color',[0.5 0.5 0.5],'units','points','extent',[mx{j}/mx{3}-0.1,(s-i+1)/s 0.1 0.1])
     case 'b_'
      a(1:2) = [];
      text(mx{j}/mx{3}-0.1,(s-i+1)/s,a,'color','b','units','normalized')
     otherwise
      a(1:2) = [];
      text(mx{j}/mx{3}-0.1,(s-i+1)/s,a,'color','m','units','normalized')
    end
  end
end
%set(gca,'xlim',[0 mx{3}+5000],'ylim',[0 s+10])
