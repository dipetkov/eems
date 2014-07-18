

function ll = wishart_ll(Sstruct,kernel,params)
%% Compute the log likelihood for               %%
%% Sv ~ Wishart(-L*Delta*L'*(s2loc/df),df)      %%
%% where Sv is the observed similarity matrix   %%
%% from microsats with possibly missing alleles %%
%% Note that because of the loop, it would be   %%
%% much faster to use sats/wishart_ll.m if      %%
%% there are no missing alleles                 %%


n = Sstruct.nIndiv;
s2loc = params.s2loc;
nSites = Sstruct.nSites;
trDinvQxD = zeros(nSites,1);
ldetDinvQ = zeros(nSites,1);

for s = 1:nSites
  X = kernel.X{s};
  XC = kernel.XC{s};
  %% oDinvo = 1'*inv(Delta)*1                    %% where J'*1*1'*J = JtOJ
  %%        = -c'*inv(w) + trace(X*J'*1*1'*J)    %% c = diag(J*J')
  %%                                             %% w = 1 (McRae's approximation)
  %%                                             %% Therefore, c'*inv(w) = n
  %% A = trace(inv(Delta)*D)
  %%   = trace(X*J'*D*J)                         %% where Jt*D*J = JtDJ
  %% B = trace(1*1'*inv(Delta)*D*inv(Delta))
  %%   = inv(w)'*J'*D*J*inv(w)                   %% where inv(w)'*J'*D*J*inv(w) = sum(sum(D))
  %%   - trace(X*(c*inv(w)'*J'*D*J + J'*D*J*inv(w)*c'))
  %%   + 1'*C*X*J'*D*J*X*C*1                     %% where c*inv(w)'*J'*D*J = J'*1*1'*D*J
  oDinvo = Sstruct.oDinvoconst(s) ...
         + sum(sum(X.*Sstruct.JtOJ{s}));
  A = sum(sum(X.*Sstruct.JtDJ{s}));
  B = Sstruct.Bconst(s) ...
    - sum(sum(X.*Sstruct.JtDJvct{s})) ...
    + sum(sum(XC'*Sstruct.JtDJ{s}*XC));
  %% ldDinvconst = log(1'*1) - logdet(J*J') + logdet(Winv)
  %%     ldCiBwi = logdet(inv(C)-B*inv(W))
  %%     where C = diag(c), W = diag(J*w) and Winv = diag(J*winv)
  ldetDinvQ(s) = Sstruct.ldDinvconst(s) ...      %% logDet(-inv(Delta)*Q)
               - kernel.ldCiBwi(s) ...
               - log(abs(oDinvo));
  trDinvQxD(s) = A - B/oDinvo;                   %% trace(inv(Delta)*Q*D)
end

%% If the degrees of freedom are not updated,   %%
%% then the log likelihood across p independent %%
%% microsates is proportional to                %%
ll = sum(ldetDinvQ) ...
   - sum((n-1).*log(s2loc)) ...
   - sum(trDinvQxD./s2loc);
ll = ll/2;
