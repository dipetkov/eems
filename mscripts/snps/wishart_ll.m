

function ll = wishart_ll(Data,kernel,params)
%% Compute the log likelihood for                 %%
%% -L*Dobs*L' ~ Wishart(-L*Dhat*L'*(s2loc/df),df) %%
%%   where Dobs is the observed distance matrix   %%
%%   from SNP data (either haploid or diploid)    %%
%%     and Dhat is the model distance matrix      %%


n = Data.nIndiv;
X = kernel.X;
XC = kernel.XC;
df = params.df;
s2loc = params.s2loc;

%% oDinvo = 1'*inv(Dhat)*1
%%        = -c'*inv(w) + trace(X*J'*1*1'*J)  %% where c = diag(J*J'),
%%                                           %%   J'*1*1'*J = JtOJ
%% A = trace(inv(Dhat)*Dobs)
%%   = trace(X*J'*Dobs*J)                    %% where J'*Dobs*J = JtDJ
%% B = trace(1*1'*inv(Dhat)*Dobs*inv(Dhat))
%%   = inv(w)'*J'*Dobs*J*inv(w)
%%   - trace(X*(c*inv(w)'*J'*Dobs*J + J'*Dobs*J*inv(w)*c'))
%%   + 1'*C*X*J'*Dobs*J*X*C*1      
%%   = inv(w)'*J'*Dobs*J*inv(w) - trace(X*cvtJtDJvct) + 1'*XC'*JtDJ*XC
oDinvo = kernel.oDinvoconst ...
       + sum(sum(X.*Data.JtOJ));
A = sum(sum(X.*Data.JtDJ));
B = kernel.Bconst ...
  - sum(sum(X.*kernel.cvtJtDJvct)) ...
  + sum(sum(XC'*Data.JtDJ*XC));
ldetDinvQ = kernel.ldDinvconst ...
          + kernel.ldCiBwi ...
          - log(abs(oDinvo));
trDinvQxD = A - B/oDinvo;
lDf2s2loc = log(df/2)-log(s2loc);

ll = df*ldetDinvQ ...
   + df*(n-1)*lDf2s2loc ...
   - df*trDinvQxD/s2loc ...
   - (df-n)*Data.ldDviQ;
ll = ll/2 ...
   - mvgammaln(df/2,n-1);

if Data.Testing
  L = Data.L;
  %% Dobs and Dhat are n-by-n matrices where n is the number of
  %% samples; X is an o-by-o matrix where o is the number of 
  %% observed demes (those with at least one sample)
  Dobs = Data.Diffs;
  Dhat = kernel.Delta;
  %% The term n*logdet(L*L')/2 is a constant that does not depends on
  %% the parameters, so is not included in ll
  ll = ll - n*logdet(L*L')/2;
  ll0 = wishpdfln(-L*Dobs*L',-L*Dhat*L'*s2loc/df,df);
  %% Check the numerical precision
  fprintf(2,'|ll0-ll|/|ll0| = %.8f\n',abs(ll0-ll)/abs(ll0));
end
