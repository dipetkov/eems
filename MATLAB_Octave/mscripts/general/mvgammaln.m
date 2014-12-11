

function x = mvgammaln(a,p)
% Compute the logarithm of the multivariate gamma function


if ~isscalar(p)
  error('mvngammaln','The second argument of mvngamma must be an integer.')
end

x = log(pi)*p*(p-1)/4 + sum(gammaln(a-(0:(p-1))/2));
