

function [proposal,pi1_pi0] = ...
    birthdeath_mVoronoi(kernel,params,mVoronoi,Sstruct,Mstruct)

%%%%%%%%%%
type = 5;%
%%%%%%%%%%

mtile0 = mVoronoi.mtiles;
negBiSize = params.negBiSize;
negBiProb = params.negBiProb;
u = rand(1);

if (mtile0==1)
  %% Since there must be at least one tile, rule out a death proposal %%
  u = 0;
end

if (u<1/2)
  %% Propose birth %%
  subtype = 1;
  mtile1 = mtile0+1;
  effect = init_effects(1,params.mEffctHalfInterval,params.mrateS2);
  center = runif_habitat(1,mVoronoi.habitat);
  mSeeds = [mVoronoi.mSeeds;center];
  mEffcts = [mVoronoi.mEffcts;effect];
  pBirth = 1/2; pDeath = 1/2;
  if (mtile0==1) pBirth = 1; end
  pi1_pi0 = log(pDeath/pBirth) + log((mtile0+negBiSize)/(mtile1/negBiProb));
else
  %% Propose death %%
  subtype = 2;
  mtile1 = mtile0-1;
  i = ceil(mtile0*rand(1));
  mSeeds = mVoronoi.mSeeds;
  mEffcts = mVoronoi.mEffcts;
  mSeeds(i,:) = [];
  mEffcts(i,:) = [];
  pBirth = 1/2; pDeath = 1/2;
  if (mtile0==2) pBirth = 1; end
  pi1_pi0 = log(pBirth/pDeath) + log((mtile0/negBiProb)/(mtile1+negBiSize));
end

proposal = struct('type',{type},'subtype',{subtype},...
                  'mEffcts',{mEffcts},...
                  'mSeeds',{mSeeds},...
                  'mtiles',{mtile1});
 
mRates = realpow(10,params.mrateMu + mEffcts);
mValues = average_rates(Mstruct,mRates,mSeeds,mVoronoi.Demes);
proposal.kernel = resistance_kernel(Sstruct,Mstruct,mValues);
