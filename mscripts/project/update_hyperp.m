

function params = update_hyperp(Sstruct,Voronoi,params,mcmc)
%% Update hyperparameters (the error variance rateS2) %%


mtiles = Voronoi.mtiles;
mEffcts = Voronoi.mEffcts;
SSm = mEffcts'*mEffcts;
a = params.mrateShape + mtiles;
b = params.mrateScale + SSm;
mrateS2 = rinvgam(a/2,b/2);
params.mrateS2 = mrateS2;
