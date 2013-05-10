function [x,firstFlag]=runOneGroupOf1000post(firstFlag,normalizedBac,fracRelevantReads,sumOneFlag)

numVariables = size(normalizedBac,2);

if sumOneFlag
  cvx_begin
    variable x(numVariables)
    minimize( norm(normalizedBac*x-fracRelevantReads) );
    subject to 
    x >= 0;
    sum( x )==1;
  cvx_end
else
  cvx_begin
    variable x(numVariables)
    minimize( norm(normalizedBac*x-fracRelevantReads) );
    subject to 
        x >= 0;
  cvx_end
end









