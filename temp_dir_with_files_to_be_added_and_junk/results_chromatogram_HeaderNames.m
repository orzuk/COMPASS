function zz=results_chromatogram_HeaderNames(zz,maxLength)
%keyboard
for i=1:length(zz)
    if iscell(zz{i,1})
        st='';
        if length(zz{i,1})>10 & length(zz{i,1})<20
            for j=1:length(zz{i,1})
                st=[st,';',zz{i,1}{j}(1:min([50,length(zz{i,1}{j})]))];
            end
        elseif length(zz{i,1})>=20
            for j=1:length(zz{i,1})
                st=[st,';',zz{i,1}{j}(1:min([30,length(zz{i,1}{j})]))];
            end
        else
            for j=1:length(zz{i,1})
                st=[st,';',zz{i,1}{j}];
            end
        end
    end
    if length(st)<maxLength
        zz{i,1}=st(2:end);
    else
        zz{i,1}=st(2:maxLength);
    end
end
