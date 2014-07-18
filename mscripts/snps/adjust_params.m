

function params = adjust_params(Sstruct,kernel,params)
%% Initialize the scale parameters (s2loc) %%


%%%%%%%%%%%%%%%%%%%%%%
params.numthetas = 2;%
%%%%%%%%%%%%%%%%%%%%%%
df = params.df;
XC = kernel.XC;
X = kernel.X;
n = Sstruct.nIndiv;
oDinvo = Sstruct.oDinvoconst ...
       + sum(sum(X.*Sstruct.JtOJ));
A = sum(sum(X.*Sstruct.JtDJ));
B = Sstruct.Bconst ...
  - sum(sum(X.*Sstruct.JtDJvct)) ...
  + sum(sum(XC'*Sstruct.JtDJ*XC));
trDinvQxD = A - B/oDinvo;
c = params.s2locShape + df*(n-1);
d = params.s2locScale + df*trDinvQxD;
s2loc = rinvgam(c/2,d/2);
params.s2loc = s2loc;
