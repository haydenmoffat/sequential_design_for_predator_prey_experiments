% This script determines the design point from the discrete design space
% which maximises the total entropy utility

Nt = Nmin:Nmax; % define possible design points
utility = zeros(length(Nt),1); % set up vector of utility values

for j = 1:length(Nt)
    clearvars log_wsum B b log_p_y log_Z_n
    
    log_wsum = zeros(K, Nt(j)+1);
    B = zeros(K, Nt(j)+1);
    
    for M = 1:K
        clearvars parameters y llh log_w_hat log_W_hat
        
        % COMPUTE THE LOG LIKELIHOOD FOR EACH POSSIBLE Y VALUE USING THE SET OF PARTICLES FROM MODEL M
        parameters = find_parameters(theta(:,:,M), Nt(j), time, Models(M));
        y = 0:Nt(j); % All possible y values
        llh = log_lik(y, Nt(j), parameters(:,1), parameters(:,2)); % Determine the log likelihood of observing each value of y at design point Nt(j) for each parameter set
        
        log_w_hat = log(W(:,M))+ llh; % updated log unnormalised weights of particles if y was the new observation
        log_wsum(M,:) = logsumexp(log_w_hat,1);
        log_W_hat = log_w_hat - log_wsum(M,:); % updated log normalised weights
        
        b = llh.*exp(log_W_hat); % multiply log likelihood by normalised weights
        b(isnan(b)) = 0;
        B(M,:) = sum(b,1);
    end
    
    log_Z_n = log_Z - logsumexp(log_Z); % normalised evidence
    log_p_y = logsumexp(log_Z_n + log_wsum.',2); % posterior predictive probabilities
    log_p_y = log_p_y - logsumexp(log_p_y); % normalised posterior predictive probabilities
    
    
    % DETERMINE THE EXPECTED UTILITY
    utility(j) = exp(log_Z_n)* sum(exp(log_wsum).* B,2) - (log_p_y.' * exp(log_p_y));
end

% PLOT THE EXPECTED UTILITIES FOR EACH DESIGN POINT
figure;
plot(Nt, utility, '.', 'Markersize', 12);
ylabel('U')
xlabel('Initial number of prey (N)')
drawnow;

% Show optimal design point
idx = find(max(utility)==utility);
disp(['Optimal design point at ', num2str(Nt(idx))]);