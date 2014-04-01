

function x = rinvgam(a,b,c)
%% Parametrize in terms of shape and scale %%
%% rather than shape and rate              %%


if (nargin==2)
  shape = a;
  scale = b;
  if isscalar(a)
    n = nrow(b);
    m = ncol(b);
  elseif isscalar(b)
    n = nrow(a);
    m = ncol(a);
  elseif ~isequal(size(a),size(b))
    error('rinvgam:WrongUsage','usage:rinvgam(shape,scale).')
  end
  x = scale./randg(shape,n,m);
elseif (nargin==3)
  if ~(isscalar(a) && isscalar(b) && isscalar(c))
    error('rinvgam:WrongUsage','usage:rinvgam(n,shape,scale).')
  end
  shape = b;
  scale = c;
  x = scale./randg(shape,a,1);
else
  error('rinvgam:NotEnoughInputs','Wrong number of parameters.')
end
