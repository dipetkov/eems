

function params = adjust_params(Sstruct,kernel,params)
%% Initialize the scale parameters (s2loc) %%


%%%%%%%%%%%%%%%%%%%%%%
params.numthetas = 1;%
%%%%%%%%%%%%%%%%%%%%%%
nSites = Sstruct.nSites;
s2loc = zeros(nSites,1);

for s = 1:nSites
  X = kernel.X{s};
  XC = kernel.XC{s};
  n = Sstruct.nIndiv(s);
  oDinvo = Sstruct.oDinvoconst(s) ...
         + sum(sum(X.*Sstruct.JtOJ{s}));
  A = sum(sum(X.*Sstruct.JtDJ{s}));
  B = Sstruct.Bconst(s) ...
    - sum(sum(X.*Sstruct.JtDOJ{s})) ...
    + sum(sum(XC'*Sstruct.JtDJ{s}*XC));
  trDinvQxD = (A - B/oDinvo)/2;
  c = params.s2locC + (n-1);
  d = params.s2locD + trDinvQxD;
  s2loc(s) = rinvgam(c/2,d/2);
end

params.s2loc = s2loc;
