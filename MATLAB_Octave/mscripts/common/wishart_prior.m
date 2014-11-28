

function logpi = wishart_prior(Data,qVoronoi,mVoronoi,params)
%% Prior probability on the log scale %%


df = params.df;
s2loc = params.s2loc;

negBiSize = params.negBiSize;
negBiProb = params.negBiProb;
mrateShape = params.mrateShape;
mrateScale = params.mrateScale;
qrateShape = params.qrateShape;
qrateScale = params.qrateScale;
s2locShape = params.s2locShape;
s2locScale = params.s2locScale;

mtiles = mVoronoi.mtiles;
qtiles = qVoronoi.qtiles;
mEffcts = mVoronoi.mEffcts;
qEffcts = qVoronoi.qEffcts;
mrateMu = params.mrateMu;
mrateS2 = params.mrateS2*2;
qrateS2 = params.qrateS2*2;

inrange = min(is_in_habitat(mVoronoi.habitat,mVoronoi.mSeeds));
inrange = (inrange && min(is_in_habitat(qVoronoi.habitat,qVoronoi.qSeeds)));
inrange = (inrange && min(abs(mEffcts)<params.mEffctHalfInterval));
inrange = (inrange && min(abs(qEffcts)<params.qEffctHalfInterval));
inrange = (inrange && abs(mrateMu)<params.mrateMuHalfInterval);
inrange = (inrange && df>=params.dfmin && df<=params.dfmax);

%% prior = NegativeBinomial(T;negBiSize,negBiProb)
%% * InverseGamma(rateS2;ratesShape,ratesScale)
%% * prod_{s=1}^S InverseGamma(s2loc_s;s2locShape,s2locScale)) %% for each microsat
%% * prod_{t=1}^T Normal(Effect_t;0,rateS2)                    %% for each Voronoi tile

if (inrange)
  logpi = - log(df) ...
	  + gammaln(negBiSize+mtiles) - gammaln(mtiles+1) + mtiles*log(negBiProb) ...
	  + gammaln(negBiSize+qtiles) - gammaln(qtiles+1) + qtiles*log(negBiProb) ...
          - ((mrateShape/2+1)*log(mrateS2) + (mrateScale/2)/mrateS2) ...
          - ((qrateShape/2+1)*log(qrateS2) + (qrateScale/2)/qrateS2) ...
          - ((mtiles/2)*log(pi*mrateS2) + (mEffcts'*mEffcts)/mrateS2) ...
          - ((qtiles/2)*log(pi*qrateS2) + (qEffcts'*qEffcts)/qrateS2) ...
	  - sum((s2locShape/2+1)*log(s2loc) + (s2locScale/2)./s2loc);
else
  logpi = -Inf;
end
