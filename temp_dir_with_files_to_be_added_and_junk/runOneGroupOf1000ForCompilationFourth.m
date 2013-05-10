function [x]=runOneGroupOf1000ForCompilationFourth(normalizedBac,fracRelevantReads)

numVariables = size(normalizedBac,2);
%keyboard
cvx_begin
  cvx_quiet(true)     
  variable x(numVariables)
  minimize( norm(normalizedBac*x-fracRelevantReads) );
  subject to 
  x >= 0;
cvx_end

x = x./sum(x);










