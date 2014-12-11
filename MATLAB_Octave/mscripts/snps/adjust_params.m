

function params = adjust_params(Data,kernel,params)
%% Initialize the scale parameter s2loc       %%
%% The other scale parameter is df, which is  %%
%% already initiailized as df = n (nIndiv)    %%


%%%%%%%%%%%%%%%%%%%%%%
params.numthetas = 2;%
%%%%%%%%%%%%%%%%%%%%%%
df = params.df;
n = Data.nIndiv;
oDinvo = kernel.oDinvo;
Xc_winv = kernel.Xc_winv;
A = sum(sum(kernel.X.*Data.JtDJ));
B = Xc_winv'*Data.JtDJ*Xc_winv;
trDinvQxD = A - B/oDinvo;
c = params.s2locShape + df*(n-1);
d = params.s2locScale + df*trDinvQxD;
s2loc = rinvgam(c/2,d/2);
params.s2loc = s2loc;
