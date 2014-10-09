

function x = rnegbin(n,r,p)
% Sample from zero truncated negative binomial with support N+ = {1,2,3,...}


if (nargin~=3)
  error('rnegbin:NotEnoughInputs','Wrong number of parameters.')
end


x = zeros(0,1);
while (length(x)<n)
  %%t = nbinrnd(r,1-p,n,1);
  t = randraw('negbinom',[r,1-p],1);
  x = [x;t];
  x = x(x>0);
end
x = x(1:n);
