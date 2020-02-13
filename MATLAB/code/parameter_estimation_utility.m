% This script calculates the expected parameter estimation utility for specified design
% points. It uses the current particle set (theta) and current normalised weights
% of the particles (W) from the SMC output

Nt = Nmin:Nmax; % possible design points

utility = zeros(length(Nt), 1);
u1 = zeros(K,1);

for j = 1:length(Nt)
    disp(j)
    for M=1:K
        clearvars parameters y llh Z_next W_hat w_hat u
        parameters = find_parameters(theta(:,:,M), Nt(j), time, Models(M));  % Determine the mu and lambda parameters for all particles
        y = 0:Nt(j); % Define possible responses
        llh = log_lik(y, Nt(j), parameters(:,1), parameters(:,2)); % Calculate the log likelihood of observing the datapoint [Nt(j), y]
        Z_next = exp(llh).'*W(:,M); % Calculate the marginal likelihood
        W_hat = exp(llh).*W(:,M); % Determine updated unnormalised weights
        w_hat = W_hat ./ sum(W_hat,1); % Normalise updated weights
        A = llh.*w_hat; % Multiply log likelihood by normalised weights
        A(isnan(A)) = 0;
        u = sum(A,1) - log(sum(W_hat,1)); % Determine utility for design point Nt(j), observation y and model M
        u1(M) = (u*Z_next); % Average utility over all responses
    end
    log_Z_n = log_Z - logsumexp(log_Z); % Normalise evidence
    utility(j) = exp(log_Z_n) * u1; % Average utility over all models
end

% PLOT THE EXPECTED UTILITIES FOR EACH DESIGN POINT
figure;
plot(Nt, utility, '.', 'Markersize', 12);
ylabel('U');
xlabel('Initial number of prey (N)');
drawnow;

% DISPLAY THE OPTIMAL DESIGN POINT
idx = find(max(utility)==utility);
if length(idx)>1
    idx = idx(1);
end
disp(['Optimal design point at ', num2str(Nt(idx))]);
