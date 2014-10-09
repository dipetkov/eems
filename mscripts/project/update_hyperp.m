

function params = update_hyperp(Data,qVoronoi,mVoronoi,params,mcmc)
%% Update hyperparameters (the effects variance rateS2) %%


mtiles = mVoronoi.mtiles;
mEffcts = mVoronoi.mEffcts;
SSm = mEffcts'*mEffcts;
a = params.mrateShape + mtiles;
b = params.mrateScale + SSm;
mrateS2 = rinvgam(a/2,b/2);
params.mrateS2 = mrateS2;

qtiles = qVoronoi.qtiles;
qEffcts = qVoronoi.qEffcts;
SSq = qEffcts'*qEffcts;
c = params.qrateShape + qtiles;
d = params.qrateScale + SSq;
qrateS2 = rinvgam(c/2,d/2);
params.qrateS2 = qrateS2;
