

function params = adjust_params(Sstruct,kernel,params)
%% Initialize the scale parameter s2loc                %%
%% This is a vector with one element for each microsat %%


%%%%%%%%%%%%%%%%%%%%%%
params.numthetas = 1;%
%%%%%%%%%%%%%%%%%%%%%%
nSites = Sstruct.nSites;
s2loc = zeros(nSites,1);

for s = 1:nSites
  X = kernel.X{s};
  XC = kernel.XC{s};
  n = Sstruct.nIndiv(s);
  oDinvo = kernel.oDinvoconst(s) ...
         + sum(sum(X.*Sstruct.JtOJ{s}));
  A = sum(sum(X.*Sstruct.JtDJ{s}));
  B = kernel.Bconst(s) ...
    - sum(sum(X.*kernel.cvtJtDJvct{s})) ...
    + sum(sum(XC'*Sstruct.JtDJ{s}*XC));
  trDinvQxD = A - B/oDinvo;
  c = params.s2locShape + (n-1);
  d = params.s2locScale + trDinvQxD;
  s2loc(s) = rinvgam(c/2,d/2);
end

params.s2loc = s2loc;
