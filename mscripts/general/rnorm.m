

function x = rnorm(Mean,Var)


if (nargin~=2)
  error('rnorm:WrongUsage','usage:rnorm(Mean,Var).')
end

if isscalar(Mean)
  sampleSize = size(Var);
elseif isscalar(Var)
  sampleSize = size(Mean);
elseif isequal(size(Mean),size(Var))
  sampleSize = size(Mean);
else
  error('rnorm:WrongUsage','usage:rnorm(Mean,Var).')
end

x = Mean+sqrt(Var).*randn(sampleSize);
