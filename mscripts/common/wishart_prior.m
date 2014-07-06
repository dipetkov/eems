

function logpi = wishart_prior(Sstruct,Voronoi,params)
%% Prior probability on the log scale  %%
%% up to a constant of proportionality %%


s2loc = params.s2loc;
negBiR = params.negBiR;
negBiP = params.negBiP;

lmlob = params.minRate;
lmupb = params.maxRate;
ratesC = params.ratesC;
ratesD = params.ratesD;
s2locC = params.s2locC;
s2locD = params.s2locD;

df = params.df;
T = Voronoi.ntiles;

Effcts = Voronoi.Effcts;
Scoord = Voronoi.Scoord;
rateMu = params.rateMu;
rateS2 = params.rateS2*2;

inrange = min(is_in_habitat(Voronoi.habitat,Scoord));
inrange = (inrange && min(abs(Effcts)<params.absDiff));
inrange = (inrange && rateMu>lmlob && rateMu<lmupb);
inrange = (inrange && df>params.dfmin && df<params.dfmax);

%% prior 
%% = NegativeBinomial(T;negBiR,negBiP)
%% * InverseGamma(rateS2;ratesC,ratesD)
%% * prod_{s=1}^S InverseGamma(s2loc_s;s2locC,s2locD)) %% for each microsatellite
%% * prod_{t=1}^T Normal(Effect_t;0,rateS2)            %% for each Voronoi tile

if (inrange)
  logpi = - log(df) ...
	  + gammaln(negBiR+T) - gammaln(T+1) + T*log(negBiP) ...
          - ((ratesC/2+1)*log(rateS2) + (ratesD/2)/rateS2) ...
          - sum((s2locC/2+1)*log(s2loc) + (s2locD/2)./s2loc) ...
          - ((T/2)*log(pi*rateS2) + (Effcts'*Effcts)/rateS2);
else
  logpi = -Inf;
end
