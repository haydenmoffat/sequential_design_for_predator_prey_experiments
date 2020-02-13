% Performs a stable calculation of log(sum(exp(x)))
function f = logsumexp(x, dim)

if (nargin < 2)
    dim = 0;    % default value for dim
end

if dim == 0
    my_max = max(x(:));
    x = x - my_max;
    f = my_max + log(sum(exp(x(:))));
else
    my_max = max(x,[],dim);
    x = x - my_max;
    f = my_max + log(sum(exp(x),dim));
end

end