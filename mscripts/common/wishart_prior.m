

function logpi = wishart_prior(Sstruct,mVoronoi,params)
%% Prior probability on the log scale  %%
%% up to a constant of proportionality %%


df = params.df;
s2loc = params.s2loc;

negBiSize = params.negBiSize;
negBiProb = params.negBiProb;
mrateShape = params.mrateShape;
mrateScale = params.mrateScale;
s2locShape = params.s2locShape;
s2locScale = params.s2locScale;

mtiles = mVoronoi.mtiles;
mSeeds = mVoronoi.mSeeds;
mEffcts = mVoronoi.mEffcts;
mrateMu = params.mrateMu;
mrateS2 = params.mrateS2*2;

inrange = min(is_in_habitat(mVoronoi.habitat,mSeeds));
inrange = (inrange && abs(mrateMu)<params.mrateMuHalfInterval);
inrange = (inrange && min(abs(mEffcts)<params.mEffctHalfInterval));
inrange = (inrange && df>params.dfmin && df<params.dfmax);

%% prior 
%% = NegativeBinomial(T;negBiSize,negBiProb)
%% * InverseGamma(mrateS2;mrateShape,mrateScale)
%% * prod_{s=1}^S InverseGamma(s2loc_s;s2locShape,s2locScale)) %% for each microsatellite
%% * prod_{t=1}^T Normal(Effect_t;0,rateS2)                    %% for each Voronoi tile

if (inrange)
  logpi = - log(df) ...
	  + gammaln(negBiSize+mtiles) - gammaln(mtiles+1) + mtiles*log(negBiProb) ...
          - ((mrateShape/2+1)*log(mrateS2) + (mrateScale/2)/mrateS2) ...
          - ((mtiles/2)*log(pi*mrateS2) + (mEffcts'*mEffcts)/mrateS2) ...
          - sum((s2locShape/2+1)*log(s2loc) + (s2locScale/2)./s2loc);
else
  logpi = -Inf;
end
