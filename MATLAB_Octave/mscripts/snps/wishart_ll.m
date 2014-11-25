

function ll = wishart_ll(Data,kernel,params)
%% Compute the log likelihood for                 %%
%% -L*Dobs*L' ~ Wishart(-L*Dhat*L'*(s2loc/df),df) %%
%%   where Dobs is the observed distance matrix   %%
%%   from SNP data (either haploid or diploid)    %%
%%     and Dhat is the model distance matrix      %%


n = Data.nIndiv;
df = params.df;
s2loc = params.s2loc;
%% oDinvo = 1'*inv(Delta)*1
%% A = trace(X*JtDJ)
%% B = (X*c-w)'*JtDJ*(X*c-w)
oDinvo = kernel.oDinvo;
Xc_winv = kernel.Xc_winv;
A = sum(sum(kernel.X.*Data.JtDJ));
B = Xc_winv'*Data.JtDJ*Xc_winv;
ldetDinvQ = kernel.ldCiBwi ...
          + log(n) - log(abs(oDinvo));
trDinvQxD = A - B/oDinvo;
lDf2s2loc = log(df/2)-log(s2loc);

ll = df*ldetDinvQ ...
   + df*(n-1)*lDf2s2loc ...
   - df*trDinvQxD/s2loc ...
   - (df-n)*Data.ldDviQ ...
   - n*Data.ldLLt;
ll = ll/2 ...
   - mvgammaln(df/2,n-1);

if Data.Testing
  L = Data.L;
  %% Dobs and Dhat are n-by-n matrices where n is the number of
  %% samples; X is an o-by-o matrix where o is the number of 
  %% observed demes (those with at least one sample)
  Dobs = Data.Diffs;
  Dhat = kernel.Delta;
  ll0 = wishpdfln(-L*Dobs*L',-L*Dhat*L'*s2loc/df,df);
  %% Check the numerical precision
  fprintf(2,'|ll0-ll|/|ll0| = %.8f\n',abs(ll0-ll)/abs(ll0));
end
