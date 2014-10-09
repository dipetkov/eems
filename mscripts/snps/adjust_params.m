

function params = adjust_params(Data,kernel,params)
%% Initialize the scale parameter s2loc       %%
%% The other scale parameter is df, which is  %%
%% already initiailized as df = n (nIndiv)    %%


%%%%%%%%%%%%%%%%%%%%%%
params.numthetas = 2;%
%%%%%%%%%%%%%%%%%%%%%%
df = params.df;
XC = kernel.XC;
X = kernel.X;
n = Data.nIndiv;
oDinvo = kernel.oDinvoconst ...
       + sum(sum(X.*Data.JtOJ));
A = sum(sum(X.*Data.JtDJ));
B = kernel.Bconst ...
  - sum(sum(X.*kernel.cvtJtDJvct)) ...
  + sum(sum(XC'*Data.JtDJ*XC));
trDinvQxD = A - B/oDinvo;
c = params.s2locShape + df*(n-1);
d = params.s2locScale + df*trDinvQxD;
s2loc = rinvgam(c/2,d/2);
params.s2loc = s2loc;
