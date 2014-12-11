

function x = rtnorm(n,Min,Max,Mean,Var)


if (nargin~=5)
  error('rtnorm:NotEnoughInputs','Wrong number of parameters.')
end
if ~(isscalar(n) && isscalar(Min) && isscalar(Max) && isscalar(Mean) && isscalar(Var))
  error('rnorm:WrongUsage','usage:rtnorm(n,Min,Max,Mean,Var).')
end

Sigma = sqrt(Var);
sampleSize = [n 1];
%%lPHI = normcdf((Min-Mean)/Sigma);
%%rPHI = normcdf((Max-Mean)/Sigma);
%%x = (lPHI+(rPHI-lPHI)*rand(a,1));
%%x = Mean + sqrt(2)*Sigma*erfinv(2*x-1);
x = randraw('normaltrunc',[Min Max Mean Sigma],sampleSize);
