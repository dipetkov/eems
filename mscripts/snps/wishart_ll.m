

function ll = wishart_ll(Sstruct,kernel,params)
%% Compute the log likelihood for                 %%
%% -L*Dobs*L' ~ Wishart(-L*Dhat*L'*(s2loc/df),df) %%
%%   where Dobs is the observed distance matrix   %%
%%   from SNP data (either haploid or diploid)    %%
%%     and Dhat is the model distance matris      %%


n = Sstruct.nIndiv;
df = params.df;
X = kernel.X;
XC = kernel.XC;
s2loc = params.s2loc;

%% oDinvo = 1'*inv(Dhat)*1                    
%%        = -c'*inv(w) + trace(X*J'*1*1'*J)  %% where c = diag(J*J'), J'*1*1'*J = JtOJ
%% A = trace(inv(Dhat)*Dobs)
%%   = trace(X*J'*Dobs*J)                    %% where J'*Dobs*J = JtDJ
%% B = trace(1*1'*inv(Dhat)*Dobs*inv(Dhat))
%%   = inv(w)'*J'*Dobs*J*inv(w)                   
%%   - trace(X*(c*inv(w)'*J'*Dobs*J + J'*Dobs*J*inv(w)*c'))
%%          %% where (c*inv(w)'*J'*Dobs*J + J'*Dobs*J*inv(w)*c') = cvtJtDJvct
%%   + 1'*C*X*J'*Dobs*J*X*C*1                     
oDinvo = kernel.oDinvoconst ...
       + sum(sum(X.*Sstruct.JtOJ));
A = sum(sum(X.*Sstruct.JtDJ));
B = kernel.Bconst ...
  - sum(sum(X.*kernel.cvtJtDJvct)) ...
  + sum(sum(XC'*Sstruct.JtDJ*XC));
%% ldDinvconst = log(1'*1) - logdet(J*J') + logdet(Winv)
%%     ldCiBwi = logdet(inv(C)-B*inv(W))
%%     where C = diag(c), W = diag(J*w) and Winv = diag(J*winv)
ldetDinvQ = kernel.ldDinvconst ...      %% logDet(-inv(Dhat)*Q)
          - kernel.ldCiBwi ...
          - log(abs(oDinvo));
trDinvQxD = A - B/oDinvo;               %% trace(inv(Dhat)*Q*Dobs)
lDf2s2loc = log(df/2)-log(s2loc);

ll = df*ldetDinvQ ...
   + df*(n-1)*lDf2s2loc ...
   - df*trDinvQxD/s2loc ...
   - (df-n)*Sstruct.ldDviQ;
ll = ll/2 ...
   - mvgammaln(df/2,n-1);
