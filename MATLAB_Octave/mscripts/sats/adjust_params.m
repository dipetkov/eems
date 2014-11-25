

function params = adjust_params(Data,kernel,params)
%% Initialize the scale parameters (s2loc) %%


%%%%%%%%%%%%%%%%%%%%%%
params.numthetas = 1;%
%%%%%%%%%%%%%%%%%%%%%%
nSites = Data.nSites;
s2loc = zeros(nSites,1);

for s = 1:nSites
  n = Data.nIndiv(s);
  oDinvo = kernel.oDinvo(s);
  Xc_winv = kernel.Xc_winv{s};
  A = sum(sum(kernel.X{s}.*Data.JtDJ{s}));
  B = Xc_winv'*Data.JtDJ{s}*Xc_winv;
  trDinvQxD = A - B/oDinvo;
  c = params.s2locShape + (n-1);
  d = params.s2locScale + trDinvQxD;
  s2loc(s) = rinvgam(c/2,d/2);
end

params.s2loc = s2loc;
