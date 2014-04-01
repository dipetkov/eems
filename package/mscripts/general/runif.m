

function x = runif(a,b,c)


if (nargin==1)
  if ~isscalar(a)
    error('runif:WrongUsage','usage:runif(n).')
  end
  x = rand(a,1);
elseif (nargin==2)
  Min = a;
  Max = b;
  if isequal(size(a),size(b))
    [n,m] = size(a);
  else
    error('runif:WrongUsage','usage:runif(Min,Max).')
  end
  x = Min+(Max-Min).*rand(n,m);
elseif (nargin==3)
  Min = b;
  Max = c;
  if ~(isscalar(a) && isscalar(b) && isscalar(c))
    error('runif:WrongUsage','usage:runif(n,Min,Max).')
  end
  x = Min+(Max-Min)*rand(a,1);
else
  error('runif:NotEnoughInputs','Wrong number of parameters.')
end
