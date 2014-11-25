

function ll = wishpdfln( X, Sigma, df)
%% Probability density function of the  %%
%% Wishart distribution with paramaters %%
%% Sigma: scale matrix                  %%
%% df: degrees of freedom               %%


n = nrow(X);
U = chol(Sigma);
ldX = logdet(X);
ldS = logdet(U,'inv');
SiX = AinvB(U,X,'inv');

ll = - (n+1)*ldX ...
     - trace(SiX) ...
     - df*(ldS-ldX) ... 
     - df*n*log(2);
ll = ll/2 ... 
     - mvgammaln(df/2,n);
