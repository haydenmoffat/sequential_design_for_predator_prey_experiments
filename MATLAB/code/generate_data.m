% This script uses a selected design point to generate a random
% observation from either the beta-binomial distribution or the binomial
% distribution depending on the model of interest

% Set new design point
d = data(i,1);

if ismember(M_true, [1 2 3 4])
    % Beta-binomial distribution
    v = ode_solver([a Th], d, time, M_true).'; % solve ode
    
    % Determine parameters
    mu = (d - v)./d;
    alpha = mu./lambda;
    beta = (1 - mu)./lambda;
    beta(beta<0) = 0;
    
    % Generate data
    data(i,2)= binornd(d, betarnd(alpha, beta));
    
elseif ismember(M_true, [5 6 7 8])
    % Binomial distribution
    v = ode_solver([a Th], d, time, M_true).'; % solve ode
    mu = (d - v)./d; % determine parameters
    data(i,2)= binornd(d, mu); % generate data
    
end
