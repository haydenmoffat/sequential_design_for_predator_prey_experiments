% FUNCTION FOR THE LOGARITHM OF THE BINOMIAL PDF
% N is the number of trials
% n is the number of successes
% p is the probability of success

function y = log_binpdf(n,N,p)
y = (gammaln(N + 1)-gammaln(n + 1)-gammaln(N - n + 1) + (n.*log(p)) + ((N-n) .* log(1-p)));
y(isnan(y)) = -Inf;
end