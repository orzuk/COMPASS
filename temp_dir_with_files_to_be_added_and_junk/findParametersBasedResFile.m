function [Nbac_in_mixture,readlen,npower,bac_dist_flag,Nread]=findParametersBasedResFile(fileName)
%keyboard
Nbac_in_mixture = fn(fileName,'mixture');% for Nbac_in_mixture
readlen = fn(fileName,'readlen');
Nread = fn(fileName,'Nreads');
if ischar(Nread) && ~isempty(findstr(Nread,'Inf'))
  Nread = Inf;
end


npower = fn(fileName,'npower');npower = npower/10;
bac_dist_flag = fn(fileName,'flag'); % for bac_dist_flag
 
if npower==0.5
  npower = '05';
elseif npower==0
  npower = '0';
end



function strout=fn(fileName,str)
a = findstr(fileName,str);
currName = fileName(a:end);
b = find(currName=='_');
if length(b)>=2
 strout = str2num(currName(b(1)+1:b(2)-1)); 
else
  strout = str2num(currName(b(1)+1:end));
end


