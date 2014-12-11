

function ans = isdistmat(x)


isdist = 1;
eps = 1e-12;
%% Check that x is nonnegative and nonzero
m = min(min(x));
if (m<0)
  isdist = 0;
else
  %% Check that the diagonal has only zeros
  d = diag(x);
  a = min(min(d));
  b = max(max(d));
  if (a~=0)||(b~=0)
    isdist = 0;
  else
    %% Check that x is conditionally negative definite
    %% That is, x has exactly one positive eigenvalues
    e = eig(x);
    negative = sum(e<-eps);
    zero = sum((e>-eps)&(e<eps));
    positive = sum(e>eps);
    if (positive~=1)||(zero~=0)
      isdist = 0;
    end
  end
end

ans = (isdist==1);
