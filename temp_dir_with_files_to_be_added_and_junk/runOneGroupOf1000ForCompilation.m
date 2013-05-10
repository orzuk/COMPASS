function [x]=runOneGroupOf1000ForCompilation(normalizedBac,fracRelevantReads,sumOneFlag)

numVariables = size(normalizedBac,2);
%keyboard
if sumOneFlag
  cvx_begin
    cvx_quiet(false)     
    variable x(numVariables)
    minimize( norm(normalizedBac*x-fracRelevantReads) );
    subject to 
    x >= 0;
    sum( x )==1;
  cvx_end
else
  cvx_begin
    cvx_quiet(false)     
    variable x(numVariables)
    minimize( norm(normalizedBac*x-fracRelevantReads) );
    subject to 
        x >= 0;
  cvx_end
end









