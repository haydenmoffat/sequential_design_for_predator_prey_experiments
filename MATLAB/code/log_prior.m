% FUNCTION FOR LOG PRIOR
% This function computes the log prior probability of a set of parameters.
% Here we consider a log normal prior on each of the parameters where
% the mean and standard deviation of the natural logarithm of the variables
% are -1.4 and 1.35, respectively.
% x is the transformed parameters.
% M is the model of interest.

function p = log_prior(x, M)
if ismember(M, [1 2])
    p = log(normpdf(x(1), -1.4, 1.35))+log(normpdf(x(2), -1.4, 1.35))+log(normpdf(x(3), -1.4, 1.35));
elseif ismember(M, [3 4])
    p = log(normpdf(x(1), -1.4, 1.35))+log(normpdf(x(3), -1.4, 1.35));
elseif ismember(M, [5 6])
    p = log(normpdf(x(1), -1.4, 1.35))+log(normpdf(x(2), -1.4, 1.35));
elseif ismember(M, [7 8])
    p = log(normpdf(x(1), -1.4, 1.35));
end