% FUNCTION FOR LOG LIKELIHOOD
% determines the log likelihood of observing the first t data points for
% for the transformed parameters x
% x is the log transformation of the variables a, Th and lambda
% data1 is the data with design points in the first column and observations
% in the second
% time is the time that the predator has access to the prey
% M is the model of interest

function lik = log_likelihood(x, data1, time, M)

x_original = exp(x); % transform the particles back to original

[N,~,ic] = unique(data1(:,1));
v = ode_solver(x_original, N, time, M); % solve the ode for each unique initial condition
V = v(:,ic);
mu_ = (data1(:,1) - V')./data1(:,1); % compute proportion of prey eaten for each initial condition

if ismember(M, 1:4)
    % Beta-Binomial Distribution
    % Determine parameters of beta binomial distribution
    alpha_ = mu_./x_original(:,3).';
    beta_ = (1 - mu_)./x_original(:,3).';
    lik = (sum(log_betabinpdf(data1(:,2),data1(:,1),alpha_,beta_),1));
    
else
    % Binomial Distribution
    lik = (sum(log_binpdf(data1(:,2),data1(:,1),mu_),1));
end
end

