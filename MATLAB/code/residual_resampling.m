% This function conducts residual resampling for when ESS becomes too low
% W is the current weights
% N is the number of particles
% theta is the current particle values
% outputs the vector th which is a vector of the resampled theta values

function th = residual_resampling(W,N,theta)
R = floor(W*N);
sumR = sum(R);
idx= repelem(1:N, R);
th(1:sumR,:) = theta(idx,:);
wt = W*N - R;
th(sumR+1:N, :) = multinomial_resampling((wt)/(sum(wt)),N,N-sumR,theta);
end