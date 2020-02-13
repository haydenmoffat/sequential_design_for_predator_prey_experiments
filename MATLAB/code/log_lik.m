% Given the expected proportion of prey eaten (mu) and over-dispersion
% parameter (lambda) for each particle, this function determines the
% log likelihood of observing data [N, y]

function lik = log_lik(y, N, mu, lambda)
if all(lambda == 0)
    lik = (log_binpdf(y,N,mu));
else
    al = mu./lambda;
    be = (1-mu)./lambda;
    lik = (log_betabinpdf(y,N,al,be));
end