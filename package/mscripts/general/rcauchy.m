

function x = rcauchy(Location,Scale)


if (nargin~=2)
  error('rnorm:WrongUsage','usage:rcauchy(Location,Scale).')
end

if isscalar(Location)
  sampleSize = size(Scale);
elseif isscalar(Scale)
  sampleSize = size(Location);
elseif isequal(size(Location),size(Scale))
  sampleSize = size(Location);
else
  error('rnorm:WrongUsage','usage:rcauchy(Location,Scale).')
end

x = Location+Scale.*tan(pi*(rand(sampleSize)-0.5));
