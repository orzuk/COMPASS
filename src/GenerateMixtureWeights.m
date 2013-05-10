% Generate weights representing proportion (summing to one)
%
% n - number of items
% distribution_string - type of model
% model_parameters - parameters specifying the model
%
function weights_vec = GenerateMixtureWeights(n, statistical_model, model_parameters)

switch statistical_model
    %%%%%%%%%%%%%%%%%% Deterministic weights %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'constant' % all items have the same frequency
        weights_vec = ones(n,1)./n;
    case {'power', 'power-law'} % power law distribution
        weights_vec = 1./((1:n).^model_parameters); %setting frequencies as the power law. model_parameters is power

    %%%%%%%%%%%%%%%%%% Stochastic weights %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'uniform' % sample frequencies from a uniform distribution
        weights_vec = rand(n,1);
    case {'simplex', 'uniform-simplex'} % sample a point uniformly from the n-1 dimensional simplex
        weights_vec = diff( [ 0 sort(rand(n-1,1)) 1 ] );
    case 'exponential' % sample from standard exponential distribution 
        weights_vec = exprnd(1, n, 1);
    case 'log-normal' % sample according to the log-normal distribution
        weights_vec = exp(randn(n,1));   
end

weights_vec = weights_vec ./ sum(weights_vec); % normalize to one
P = randperm(n); weights_vec = weights_vec(P); % Permute elements

