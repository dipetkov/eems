

function params = update_dfsupp(Sstruct,params,mcmc)
%% During the burn-in phase, the degrees of freedom parameter
%% is kept fixed at its lowest possible value, nIndiv
%% (Actually, for 3/4-ths of the burn-in iterations)
%% In practice, this seems to improve convergence as it is
%% less likely that the chain gets stuck at an unfavorable
%% initial state


dfmin = params.dfmin;
dfmax = params.dfmax;
dflob = params.dflob;
dfupb = params.dfupb;
if mcmc.currIter>(mcmc.numBurnIter/3)
  prop = 1;
else 
  prop = 0;
end
dfmax = dfmax+prop*(dfupb-dfmax);
dfmin = dfmin+prop*(dflob-dfmin);
params.dfmax = dfmax;
params.dfmin = dfmin;
