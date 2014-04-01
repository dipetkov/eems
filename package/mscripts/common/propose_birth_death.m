

function [proposal,pi1_pi0] = ...
    propose_birth_death(kernel,params,Voronoi,Sstruct,Mstruct)

%%%%%%%%%%
type = 5;%
%%%%%%%%%%

T = Voronoi.ntiles;
negBiR = params.negBiR;
negBiP = params.negBiP;
u = runif(1);

if (T==1)
  %% Since there must be at least one tile, rule out a death proposal %%
  u = 0;
end

if (u<1/2)
  %% Propose birth %%
  ntiles = T+1;
  subtype = 1;
  effect = init_effects(1,params.absDiff,params.rateS2);
  center = runif_habitat(1,Voronoi.habitat);
  Scoord = [Voronoi.Scoord;center];
  Effcts = [Voronoi.Effcts;effect];
  pBirth = 1/2; pDeath = 1/2;
  if (T==1) pBirth = 1; end
  pi1_pi0 = log(pDeath/pBirth) + log((negBiP/(T+1))*(T+negBiR));
else
  %% Propose death %%
  ntiles = T-1;
  subtype = 2;
  i = unidrnd(T);
  Scoord = Voronoi.Scoord;
  Effcts = Voronoi.Effcts;
  Scoord(i,:) = [];
  Effcts(i,:) = [];
  pBirth = 1/2; pDeath = 1/2;
  if (T==2) pBirth = 1; end
  pi1_pi0 = log(pBirth/pDeath) + log((T/negBiP)/(T+negBiR-1));
end
 
mRates = realpow(10,params.rateMu + Effcts);
mValues = average_rates(Mstruct,mRates,Scoord,Voronoi.Vcoord);
proposal = resistance_kernel(Sstruct,Mstruct,mValues);
proposal.type = type;
proposal.subtype = subtype;
proposal.Scoord = Scoord;
proposal.Effcts = Effcts;
proposal.ntiles = ntiles;
