

function ll = pseudowishpdfln( X, Sigma, df)


n = nrow(X);
r = rank(X);
U = chol(Sigma);
if (r<n)
  if (r~=df)
    error('singwishpdfln:rank(X)','rank(X) != df.');
  end
end
ldS = logdet(U,'inv');
SiX = AinvB(U,X,'inv');
ldX = pseudo_logdet(X);
q = min(df,n);
ll = - (n+1)*ldX ...
     - trace(SiX) ...
     - df*(ldS-ldX) ...
     - df*n*log(2) ...
     - df*(n-q)*log(pi);
ll = ll/2 ...
     - mvgammaln(df/2,q);
