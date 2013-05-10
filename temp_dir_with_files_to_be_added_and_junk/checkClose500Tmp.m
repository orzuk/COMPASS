function [currResOrig,x]=checkClose500Tmp(normalizedBac)

x = rand(size(normalizedBac,2),1); 
%x(find(x>0.8)) = 0;
x=x./sum(x);

y = normalizedBac*x;

[currResOrig,firstFlag]=runOneGroupOf1000post(1,normalizedBac,y,1);
plot(currResOrig,x,'.')



