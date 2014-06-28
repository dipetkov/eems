

function [ll,trDinvQxD] = wishart_ll(Sstruct,kernel,params)
%% Compute the log likelihood for                %%
%% -L*Dv*L' ~ Wishart(-L*Delta*L'*(s2loc/df),df) %%
%%   where Dv is the observed distance matrix    %%
%%   from SNP data (either haploid or diploid)   %%


n = Sstruct.nIndiv;
df = params.df;
X = kernel.X;
XC = kernel.XC;
s2loc = params.s2loc;

oDinvo = Sstruct.oDinvoconst ...
       + sum(sum(X.*Sstruct.JtOJ));
A = sum(sum(X.*Sstruct.JtDJ));
B = Sstruct.Bconst ...
  - sum(sum(X.*Sstruct.JtDOJ)) ...
  + sum(sum(XC'*Sstruct.JtDJ*XC));
ldetDinvQ = Sstruct.ldDinvconst ...
          + kernel.ldDinv ...
          - log(abs(oDinvo));
trDinvQxD = (A - B/oDinvo)/2;
lDf2s2loc = log(df/2)-log(s2loc);

ll = df*ldetDinvQ ...
   + df*(n-1)*lDf2s2loc ...
   - df*trDinvQxD/s2loc ...
   - (df-n)*Sstruct.ldDviQ;
ll = ll/2 ...
   - mvgammaln(df/2,n-1);
