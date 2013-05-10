% Generate weights representing proportion (summing to one)
% 
% n - number of items 
% distribution_string - type of model 
% model_parameters - parameters specifying the model 
%
function weights_vec = GenerateMixtureWeights(n, statistical_model, model_parameters)

switch statistical_model
    case 'uniform'
        weights_vec = ones(n,1)./n; % all items have the same frequency
        
    case {'power', 'power-law'} % power law distribution 
        weights_vec = 1./((1:n).^model_parameters); %setting frequencies as the power law. model_parameters is power

    case 'exponential' 
        
        
end

weights_vec = weights_vec ./ sum(weights_vec); % normalize to one 
P = randperm(n); weights_vec = weights_vec(P); % Permute elements

