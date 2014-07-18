

function ll = wishart_ll(Sstruct,kernel,params)
%% Compute the log likelihood for                %%
%% -L*Dv*L' ~ Wishart(-L*Delta*L'*(s2loc/df),df) %%
%%   where Dv is the observed distance matrix    %%
%%   from SNP data (either haploid or diploid)   %%


n = Sstruct.nIndiv;
df = params.df;
X = kernel.X;
XC = kernel.XC;
s2loc = params.s2loc;

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
oDinvo = Sstruct.oDinvoconst ...
       + sum(sum(X.*Sstruct.JtOJ));
A = sum(sum(X.*Sstruct.JtDJ));
B = Sstruct.Bconst ...
  - sum(sum(X.*Sstruct.JtDJvct)) ...
  + sum(sum(XC'*Sstruct.JtDJ*XC));
%% ldDinvconst = log(1'*1) - logdet(J*J') + logdet(Winv)
%%     ldCiBwi = logdet(inv(C)-B*inv(W))
%%     where C = diag(c), W = diag(J*w) and Winv = diag(J*winv)
ldetDinvQ = Sstruct.ldDinvconst ...
          - kernel.ldCiBwi ...
          - log(abs(oDinvo));
trDinvQxD = A - B/oDinvo;
lDf2s2loc = log(df/2)-log(s2loc);

ll = df*ldetDinvQ ...
   + df*(n-1)*lDf2s2loc ...
   - df*trDinvQxD/s2loc ...
   - (df-n)*Sstruct.ldDviQ;
ll = ll/2 ...
   - mvgammaln(df/2,n-1);
