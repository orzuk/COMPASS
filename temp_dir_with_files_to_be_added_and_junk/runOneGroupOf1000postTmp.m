function [x,firstFlag]=runOneGroupOf1000postTmp(firstFlag,normalizedBac,fracRelevantReads,sumOneFlag)

numVariables = size(normalizedBac,2);

%keyboard

%if firstFlag & ~isdeployed
%  rmpath(genpath([userDir,'/fenoam/antiGchip/replicateTest/spider/']));
%  addpath(genpath([userDir,'/CS/BAC/cvx/']))
%  firstFlag = 0;
%end
%keyboard


if sumOneFlag
  cvx_begin
  cvx_precision([0 0 eps^4])
  cvx_slvtol = 1e-30
  cvx_slvitr = 100;
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









