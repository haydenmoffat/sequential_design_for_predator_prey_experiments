function th = multinomial_resampling(W,N,n,theta)
% This function conducts multinomial resampling for when ESS becomes too low
% W is the current weightings
% N is the number of particles
% n is the number of samples we require
% theta is the current theta values
% outputs the vector th which is a vector of the resampled theta values
idx= repelem(1:N, mnrnd(n, W));
th = theta(idx,:);
end