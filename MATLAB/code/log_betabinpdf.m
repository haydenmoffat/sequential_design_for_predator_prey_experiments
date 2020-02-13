% FUNCTION FOR THE LOGARITHM OF THE BETA BINOMIAL PDF
% n is the number of successes
% N is the number of trials
% al and be are the alpha and beta parameters of the beta-binomial
% distribution, respectively

function y = log_betabinpdf(n,N,al,be)

% Set beta values < 0 to 0 so that no errors are produced
be(be < 0) = 0;

y = (gammaln(N + 1)-gammaln(n + 1)-gammaln(N - n + 1) + betaln((al + n),(be + N - n)) - betaln(al,be+1000-1000));

% In case both betaln functions above evaluate to Inf (if al == 0 and n == 0),
% set log-likelihood to -Inf
y(isnan(y)) = -Inf;

end
