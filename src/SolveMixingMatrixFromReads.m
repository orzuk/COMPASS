% Solve L2 minimization
% ReadBySpeciesMat - Matrix A
% values - read sequences for each row of A
% uniqueReads - reads appearing in the next-gen data
% uniqueReads_counts - number of read from each unique type
% sumOneFlag - force normalization: \sum_i x_i=1 in the solution
%
% Output:
% x - solution
%
function x = SolveMixingMatrixFromReads( ...
    ReadBySpeciesMat, values, uniqueReads, uniqueReads_counts, sumOneFlag)

numVariables = size(ReadBySpeciesMat,2); % number of variables in the solution 
[~,i1,i2] = intersect(values, uniqueReads, 'rows');
numRelevantReads = zeros(size(values,1),1); % y vector for optimization
clear values
numRelevantReads(i1) = uniqueReads_counts(i2);
clear uniqueReads uniqueReads_counts % save memory

% % % if firstFlag % Add path 
% % %     rmpath(genpath('/home/csfaculty/shental/Weizmann/ph2users/fenoam/antiGchip/replicateTest/spider/'));
% % %     addpath(genpath('/home/csfaculty/shental/CS/BAC/cvx/'))
% % %     firstFlag = 0;
% % % end

%keyboard
fracRelevantReads = numRelevantReads./sum(numRelevantReads); % normlize vector y to sum to one

if sumOneFlag % normalize solution to one
    cvx_begin
    variable x(numVariables) % variable and length. The solution is stored in x
    minimize( norm(ReadBySpeciesMat*x-fracRelevantReads) ); % minimize L2 norm
    subject to
    x >= 0; % vector of probabilities
    sum( x )==1; % sum to one
    cvx_end
else % solve un-normalized 
    cvx_begin
    variable x(numVariables)
    minimize( norm(ReadBySpeciesMat*x-fracRelevantReads) );
    subject to
    x >= 0;
    cvx_end
end









