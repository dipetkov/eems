

function ll = wishart_ll(Data,kernel,params)
%% Compute the log likelihood for               %%
%% Sv ~ Wishart(-L*Delta*L'*(s2loc/df),df)      %%
%% where Sv is the observed similarity matrix   %%
%% from microsats with possibly missing alleles %%
%% Note that because of the loop, it would be   %%
%% much faster to use sats/wishart_ll.m if      %%
%% there are no missing alleles                 %%


n = Data.nIndiv;
s2loc = params.s2loc;
nSites = Data.nSites;
trDinvQxD = zeros(nSites,1);
ldetDinvQ = zeros(nSites,1);

%% oDinvo = 1'*inv(Delta)*1
%% A = trace(X*JtDJ)
%% B = (X*c-w)'*JtDJ*(X*c-w)
for s = 1:nSites
  oDinvo = kernel.oDinvo(s);
  Xc_winv = kernel.Xc_winv{s};
  A = sum(sum(kernel.X{s}.*Data.JtDJ{s}));
  B = Xc_winv'*Data.JtDJ{s}*Xc_winv;
  ldetDinvQ(s) = kernel.ldCiBwi(s) ...
               + log(n(s)) - log(abs(oDinvo));
  trDinvQxD(s) = A - B/oDinvo;
end

ll = sum(ldetDinvQ) ...
   - sum((n-1).*log(s2loc)) ...
   - sum(trDinvQxD./s2loc) ...
   - Data.llconst;
ll = ll/2;

if Data.Testing
  Delta = kernel.Delta;
  ll0 = 0;
  for s = 1:nSites
    n = Data.nIndiv(s);
    o = Data.O{s};
    z = Data.Z{s};
    L = [-ones(n-1,1),eye(n-1)];
    LDeltaLt = L*Delta(o,o)*L';
    ll0 = ll0 + pseudowishpdfln(2*L*z*z'*L',-LDeltaLt*s2loc(s),1);
  end
  %% Check the numerical precision
  fprintf(2,'|ll0-ll|/|ll0| = %.8f\n',abs(ll0-ll)/abs(ll0));
end
