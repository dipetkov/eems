

function params = adjust_params(Data,kernel,params)
%% Initialize the scale parameters (s2loc) %%


%%%%%%%%%%%%%%%%%%%%%%
params.numthetas = 1;%
%%%%%%%%%%%%%%%%%%%%%%
nSites = Data.nSites;
s2loc = zeros(nSites,1);

for s = 1:nSites
  X = kernel.X{s};
  XC = kernel.XC{s};
  n = Data.nIndiv(s);
  oDinvo = kernel.oDinvoconst(s) ...
         + sum(sum(X.*Data.JtOJ{s}));
  A = sum(sum(X.*Data.JtDJ{s}));
  B = kernel.Bconst(s) ...
    - sum(sum(X.*kernel.cvtJtDJvct{s})) ...
    + sum(sum(XC'*Data.JtDJ{s}*XC));
  trDinvQxD = A - B/oDinvo;
  c = params.s2locShape + (n-1);
  d = params.s2locScale + trDinvQxD;
  s2loc(s) = rinvgam(c/2,d/2);
end

params.s2loc = s2loc;
