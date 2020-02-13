% This script determines the design point from the discrete design space
% which maximises the model discrimination utility

Nt = Nmin:Nmax; % possible design points

utility = zeros(length(Nt),1);

for j = 1:length(Nt)
    
    clear parameters y llh log_w_new log_f W_new log_Z_new log_Z_new_n log_Z_n
    log_f = zeros(K, length(0:Nt(j)));
    log_Z_new = zeros(K, length(0:Nt(j)));
    
    for M = 1:K
        parameters = find_parameters(theta(:,:,M), Nt(j), time, Models(M));  % Determine the mu and lambda parameters for all the particles
        y = 0:Nt(j); % Define possible responses
        llh = log_lik(y, Nt(j), parameters(:,1), parameters(:,2)); % Determine the log likelihood of observing each value of y at design point Nt(j) for each parameter set
        
        log_w_new = log(W(:,M))+ llh; % Determine weights of particles if y was the new observation
        log_f(M,:) = logsumexp(log_w_new, 1); % Ratio of evidences
        log_Z_new(M,:) = log_Z(M) + log_f(M,:); % Updated evidence
        
    end
    log_Z_n = log_Z - logsumexp(log_Z); % Normalise current evidence
    log_Z_new_n = log_Z_new - logsumexp(log_Z_new, 1); % Normalise updated evidence
    utility(j) = exp(log_Z_n) *(sum(log_Z_new_n.*exp(log_f),2)); %  Calculate utility for design point Nt(j)
end

% PLOT THE EXPECTED UTILITIES FOR EACH DESIGN POINT
figure;
plot(Nt, utility, '.', 'Markersize', 12);
ylabel('U')
xlabel('Initial number of prey (N)')
drawnow;

% DISPLAY THE OPTIMAL DESIGN POINT
idx = find(max(utility)==utility);
if length(idx)>1
    idx = idx(1);
end
disp(['Optimal design point at ', num2str(Nt(idx))]);
