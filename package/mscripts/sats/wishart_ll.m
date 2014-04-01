

function [ll,trDinvQxD] = wishart_ll(Sstruct,kernel,params)
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

for s = 1:nSites;

  X = kernel.X{s};
  XC = kernel.XC{s};
  oDinvo = Sstruct.oDinvoconst(s) ...
         + sum(sum(X.*Sstruct.JtOJ{s}));
  A = sum(sum(X.*Sstruct.JtDJ{s}));
  B = Sstruct.Bconst(s) ...
    - sum(sum(X.*Sstruct.JtDOJ{s})) ...
    + sum(sum(XC'*Sstruct.JtDJ{s}*XC));
  ldetDinvQ(s) = Sstruct.ldDinvconst(s) ...
               + kernel.ldDinv(s) ...
               - log(abs(oDinvo));
  trDinvQxD(s) = (A - B/oDinvo)/2;

end

%% If the degrees of freedom are not updated,   %%
%% then the log likelihood across p independent %%
%% microsates is proportional to                %%
ll = sum(ldetDinvQ) ...
   - sum((n-1).*log(s2loc)) ...
   - sum(trDinvQxD./s2loc);
ll = ll/2;
