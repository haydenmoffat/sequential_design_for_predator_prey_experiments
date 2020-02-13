% For each particle in the given particle set (theta) and an initial condition (N),
% this function determines the parameter for the expected proportion
% of prey eaten (mu) and over-dispersion parameter (lambda) for a specified
% time (time) and model (M)

function para = find_parameters(theta, N, time, M)
theta_original = exp(theta); % transform the particle values

V = ode_solver(theta_original, N, time, M); % solve the ode

% Calculate mu and lambda
mu = (N - V)./N;
lambda = theta_original(:,3,:);
para = [mu lambda];
end

