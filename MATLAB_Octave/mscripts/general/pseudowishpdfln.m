

function ll = pseudowishpdfln( X, Sigma, df)


n = nrow(X);
r = rank(X);
q = n;
if (r<n)
  % If X is rank-deficient, df should be equal to the rank of X
  if (df~=r)
    error('pseudowishpdfln: rank(X)<nrow(X) but df!=rank(X)');
  end
  q = r;
else
  % If X is full-rank, then df should be at least equal to n-1   
  if (df<n)
    error('pseudowishpdfln: rank(X)=nrow(X) but df<nrow(X)')
  end
end
U = chol(Sigma);
ldS = logdet(U,'inv');
SiX = AinvB(U,X,'inv');
ldX = pseudologdet(X);
ll = - (n+1)*ldX ...
     - trace(SiX) ...
     - df*(ldS-ldX) ...
     - df*n*log(2) ...
     - df*(n-q)*log(pi);
ll = ll/2 ...
     - mvgammaln(df/2,q);
