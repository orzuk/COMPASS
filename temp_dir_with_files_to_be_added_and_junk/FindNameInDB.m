% first load 'bac16s_primers750' from Data dir
function [pos]=FindNameInDB(name,headers)
pos=0;
namelen=length(name);
for a=1:length(headers)
    b=headers{a};
    if (iscell(b))
        for c=1:length(b)
            d=b{c};
            if (sum(d(1:namelen)==name)==namelen)
                disp(a);
                disp(d);
                disp(c);
                pos=a;
            end
        end
    else
        d=b;
        if (sum(d(1:namelen)==name)==namelen)
            disp(a);
            disp(d);
            pos=a;
        end
    end 
end

            
            