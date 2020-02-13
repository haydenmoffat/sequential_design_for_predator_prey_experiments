% This script contains code to conduct a move step on the particles
% after the resampling step

% Compute parameters of MCMC kernel qt, i.e. weighted sample
% covariance matrix
cov_matrix = cov(theta(:,:,M));
cov_matrix(isnan(cov_matrix)) = 0;

% Move each particle theta with MCMC kernel for 1 iteration
n_moves = zeros(N,K); % variable that describes the number of particles that have moved after 1 iteration
for l= 1:N
    u=rand; x= theta(l,:, M);
    logpri(l,M) = log_prior(x, Models(M)); % determine log prior of initial theta values
    loglik(l,M) = log_likelihood(x, data_subset, time, Models(M)); % determine log likelihood of initial theta values
    px(l,M) = loglik(l,M)+logpri(l,M);
    
    xs=mvnrnd(x,cov_matrix); % new sample xs based on existing x from proposal pdf.
    
    % Determine log prior and log likelihood of proposed theta values
    logprior_xs = log_prior(xs, Models(M));
    loglik_xs = log_likelihood(xs, data_subset, time, Models(M));
    pxs=loglik_xs +logprior_xs;
    
    if u<min(1,exp(pxs-px(l,M)))
        theta(l,:,M)=xs;
        logpri(l,M)= logprior_xs;
        loglik(l,M)= loglik_xs;
        px(l,M) = pxs;
        n_moves(l,M) = n_moves(l,M)+1;
    end
    
end

% Determine acceptance probability and an appropriate number of iterations
% of the MCMC
prob = sum(n_moves(:,M))/N;
Rt = ceil(log(0.01)/log(1-prob));


% Displays the acceptance probability and the number of iterations
% (Rt) for each move step
fprintf('At i= %i, Acceptance probability = %.3f, Rt = %i \n', i, prob, Rt)

% Move particles the rest of the iterations
for l= 1:N
    for q=2:Rt
        u=rand; x= theta(l,:,M);
        xs=mvnrnd(x,cov_matrix); % new sample xs based on existing x from proposal pdf.
        % determine log prior and log likelihood of proposed theta values
        logprior_xs = log_prior(xs, Models(M)); 
        loglik_xs = log_likelihood(xs, data_subset, time, Models(M));
        pxs=loglik_xs +logprior_xs;
        
        if u<min(1,exp(pxs-px(l)))
            theta(l,:,M)=xs;
            logpri(l,M)= logprior_xs;
            loglik(l,M)= loglik_xs;
            px(l,M) = pxs;
            
        end
        
    end
end