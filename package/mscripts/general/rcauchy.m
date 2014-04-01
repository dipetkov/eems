

function x = rcauchy(a,b,c)


if (nargin==3)
  if ~(isscalar(a) && isscalar(b) && isscalar(c))
    error('rcauchy:WrongUsage','usage:rcauchy(n,Location,Scale).')
  end
  sampleSize = [a,1];
  Location = b;
  Scale = c;
  x = Location+Scale*tan(pi*(rand(sampleSize)-0.5));
elseif (nargin==2)
  Location = a;
  Scale = b;
  if isscalar(a)
    sampleSize = size(b);
  elseif isscalar(b)
    sampleSize = size(a);
  elseif isequal(size(a),size(b))
    sampleSize = size(a);
  else
    error('rcauchy:WrongUsage','usage:rcauchy(Location,Scale).')
  end
  x = Location+Scale.*tan(pi*(rand(sampleSize)-0.5));
elseif (nargin==1)
  if isscalar(a)
    sampleSize = [a,1];
    x = tan(pi*(rand(sampleSize)-0.5));
  else
    error('rnorm:WrongUsage','usage:rcauchy(n).')
  end
else
  error('rnorm:NotEnoughInputs','Wrong number of parameters.')
end
