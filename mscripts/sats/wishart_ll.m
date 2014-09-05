

function ll = wishart_ll(Sstruct,kernel,params)
%% Compute the log likelihood for                 %%
%% -L*Dobs*L' ~ Wishart(-L*Dhat*L'*(s2loc/df),df) %%
%% where Sv is the observed similarity matrix     %%
%% from microsats with possibly missing alleles   %%
%% This function loops over the microsats         %%


n = Sstruct.nIndiv;
s2loc = params.s2loc;
nSites = Sstruct.nSites;
trDinvQxD = zeros(nSites,1);
ldetDinvQ = zeros(nSites,1);

%% oDinvo = 1'*inv(Dhat)*1                    
%%        = -c'*inv(w) + trace(X*J'*1*1'*J)  %% where c = diag(J*J'), J'*1*1'*J = JtOJ
%% A = trace(inv(Dhat)*Dobs)
%%   = trace(X*J'*Dobs*J)                    %% where J'*Dobs*J = JtDJ
%% B = trace(1*1'*inv(Dhat)*Dobs*inv(Dhat))
%%   = inv(w)'*J'*Dobs*J*inv(w)                   
%%   - trace(X*(c*inv(w)'*J'*Dobs*J + J'*Dobs*J*inv(w)*c'))
%%          %% where (c*inv(w)'*J'*Dobs*J + J'*Dobs*J*inv(w)*c') = cvtJtDJvct
%%   + 1'*C*X*J'*Dobs*J*X*C*1                     
for s = 1:nSites
  X = kernel.X{s};
  XC = kernel.XC{s};
  oDinvo = kernel.oDinvoconst(s) ...
         + sum(sum(X.*Sstruct.JtOJ{s}));
  A = sum(sum(X.*Sstruct.JtDJ{s}));
  B = kernel.Bconst(s) ...
    - sum(sum(X.*kernel.cvtJtDJvct{s})) ...
    + sum(sum(XC'*Sstruct.JtDJ{s}*XC));
  %% ldDinvconst = log(1'*1) - logdet(J*J') + logdet(Winv)
  %%     ldCiBwi = logdet(inv(C)-B*inv(W))
  %%     where C = diag(c), W = diag(J*w) and Winv = diag(J*winv)
  ldetDinvQ(s) = kernel.ldDinvconst(s) ...   %% logDet(-inv(Dhat)*Q)
               - kernel.ldCiBwi(s) ...
               - log(abs(oDinvo));
  trDinvQxD(s) = A - B/oDinvo;               %% trace(inv(Dhat)*Q*Dobs)
end

%% If the degrees of freedom are not updated,    %%
%% then the log likelihood across p independent  %%
%% microsates is proportional to                 %%
ll = sum(ldetDinvQ) ...
   - sum((n-1).*log(s2loc)) ...
   - sum(trDinvQxD./s2loc);
ll = ll/2;
