

function x = rtnorm(a,b,c,d,e)


if (nargin==4)
  if ~(isscalar(a)&&isscalar(b)) && ~isequal(size(a),size(b),size(c),size(d))
    error('rtnorm:WrongUsage','usage:rtnorm(Min,Max,Mean,Var).')
  end
  Min = a;
  Max = b;
  Mean = c;
  Sigma = sqrt(d);
  if isscalar(c)
    n = nrow(d);
    m = ncol(d);
  elseif isscalar(d)
    n = nrow(c);
    m = ncol(c);
  elseif isequal(size(c),size(d))
    [n,m] = size(c);
  else
    error('rtnorm:WrongUsage','usage:rtnorm(n,Min,Max,Mean,Var).')
  end
  lPHI = normcdf((Min-Mean)./Sigma);
  rPHI = normcdf((Max-Mean)./Sigma);
  x = (lPHI+(rPHI-lPHI).*rand(n,m));
  x = Mean + sqrt(2)*Sigma.*erfinv(2*x-1);
elseif (nargin==5)
  Min = b;
  Max = c;
  Mean = d;
  Sigma = sqrt(e);
  if ~(isscalar(a) && isscalar(b) && isscalar(c) && isscalar(d) && isscalar(e))
    error('rnorm:WrongUsage','usage:rnorm(n,Mean,Var).')
  end
  lPHI = normcdf((Min-Mean)/Sigma);
  rPHI = normcdf((Max-Mean)/Sigma);
  x = (lPHI+(rPHI-lPHI)*rand(a,1));
  x = Mean + sqrt(2)*Sigma*erfinv(2*x-1);
elseif (nargin~=5)
  error('rtnorm:NotEnoughInputs','Wrong number of parameters.')
end

