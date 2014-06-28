

function params = update_hyperp(Sstruct,Voronoi,params,mcmc)
%% Update hyperparameters (the error variance rateS2) %%


ntiles = Voronoi.ntiles;
Effcts = Voronoi.Effcts;
SSe = Effcts'*Effcts;
a = params.ratesC + ntiles;
b = params.ratesD + SSe;
rateS2 = rinvgam(a/2,b/2);
params.rateS2 = rateS2;
