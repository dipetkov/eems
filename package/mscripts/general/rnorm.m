

function x = rnorm(a,b,c)


if (nargin==3)
  if ~(isscalar(a) && isscalar(b) && isscalar(c))
    error('rnorm:WrongUsage','usage:rnorm(n,Mean,Var).')
  end
  Mean = b;
  Var = c;
  x = Mean+sqrt(Var)*randn(a,1);
elseif (nargin==2)
  Mean = a;
  Var = b;
  if isscalar(a)
    n = nrow(b);
    m = ncol(b);
  elseif isscalar(b)
    n = nrow(a);
    m = ncol(a);
  elseif isequal(size(a),size(b))
    [n,m] = size(a);
  else
    error('rnorm:WrongUsage','usage:rnorm(Mean,Var).')
  end
  x = Mean+sqrt(Var).*randn(n,m);
else
  error('rnorm:NotEnoughInputs','Wrong number of parameters.')
end
