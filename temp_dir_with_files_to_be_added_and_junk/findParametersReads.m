function [Nbac_in_mixture,readlen,npower,bac_dist_flag,Nread]=findParametersReads(fileName)

Nbac_in_mixture = fn(fileName,'Nbacmix');
readlen = fn(fileName,'Readlen');
Nread = fn(fileName,'Nread');
npower = fn(fileName,'Npower');npower = npower/10;
bac_dist_flag = fn(fileName,'bacdistflag');
 

function strout=fn(fileName,str)
a = findstr(fileName,str);
currName = fileName(a:end);
b = find(currName=='_');
strout = str2num(currName(b(1)+1:b(2)-1));

