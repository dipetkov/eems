

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
for s = 1:nSites
  X = kernel.X{s};
  XC = kernel.XC{s};
  oDinvo = kernel.oDinvoconst(s) ...
         + sum(sum(X.*Data.JtOJ{s}));
  A = sum(sum(X.*Data.JtDJ{s}));
  B = kernel.Bconst(s) ...
    - sum(sum(X.*kernel.cvtJtDJvct{s})) ...
    + sum(sum(XC'*Data.JtDJ{s}*XC));
  ldetDinvQ(s) = kernel.ldDinvconst(s) ...
               - kernel.ldCiBwi(s) ...
               - log(abs(oDinvo));
  trDinvQxD(s) = A - B/oDinvo;
end

%% Since the degrees of freedom are not updated in this case, %%
%% then the log likelihood across p independent microsates is %%
%% proportional to:                                           %%
ll = sum(ldetDinvQ) ...
   - sum((n-1).*log(s2loc)) ...
   - sum(trDinvQxD./s2loc);
ll = ll/2;

if Data.Testing
  Delta = kernel.Delta;
  ll0 = 0;
  for s = 1:nSites
    n = Data.nIndiv(s);
    o = Data.O{s};
    Z = Data.Z{s};
    L = [-ones(n-1,1),eye(n-1)];
    LDeltaLt = L*Delta(o,o)*L';
    LDLt = -2*L*Z*Z'*L';
    ll0 = ll0 + pseudowishpdfln(-LDLt,-LDeltaLt*s2loc(s),1);
    ll = ll - logdet(L*L')/2 ...
       - (n-1)*pseudo_logdet(-LDLt,1)/2 ...
       - (n-2)*log(pi)/2 - (n-1)*log(2)/2 ...
       - mvgammaln(1/2,1);
  end
  %% Check the numerical precision
  fprintf(2,'|ll0-ll|/|ll0| = %.8f\n',abs(ll0-ll)/abs(ll0));
end
