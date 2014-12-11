

function [proposal,pi1_pi0] = ...
    birthdeath_qVoronoi(kernel,params,qVoronoi,mVoronoi,Data,Graph)

%%%%%%%%%%
type = 3;%
%%%%%%%%%%

qtile0 = qVoronoi.qtiles;
qSeeds = qVoronoi.qSeeds;
mSeeds = mVoronoi.mSeeds;
qEffcts = qVoronoi.qEffcts;
mEffcts = mVoronoi.mEffcts;
negBiSize = params.negBiSize;
negBiProb = params.negBiProb;

u = rand(1);

if (qtile0==1)
  %% Since there must be at least one tile, rule out a death proposal %%
  u = 0;
end

if (u<1/2)
  %% Propose birth %%
  subtype = 1;
  qtile1 = qtile0+1;
  effect = init_effects(1,params.qEffctHalfInterval,params.qrateS2);
  center = runif_habitat(1,qVoronoi.habitat);
  qSeeds = [qSeeds;center];
  qEffcts = [qEffcts;effect];
  pBirth = 1/2; pDeath = 1/2;
  if (qtile0==1) pBirth = 1; end
  pi1_pi0 = log(pDeath/pBirth) + log((qtile0+negBiSize)/(qtile1/negBiProb));
else
  %% Propose death %%
  subtype = 2;
  qtile1 = qtile0-1;
  i = ceil(qtile0*rand(1));
  qSeeds(i,:) = [];
  qEffcts(i,:) = [];
  pBirth = 1/2; pDeath = 1/2;
  if (qtile0==2) pBirth = 1; end
  pi1_pi0 = log(pBirth/pDeath) + log((qtile0/negBiProb)/(qtile1+negBiSize));
end

proposal = struct('type',{type},'subtype',{subtype},...
                  'qEffcts',{qEffcts},...
                  'qSeeds',{qSeeds},...
                  'qtiles',{qtile1});

mRates = realpow(10,mEffcts + params.mrateMu);
qRates = realpow(10,qEffcts); % qrateMu = 0.0;
[qValues,mValues] = ...
  average_rates(Graph,qRates,mRates,qSeeds,mSeeds,qVoronoi.Demes);
proposal.kernel = resistance_kernel(Data,Graph,mValues,qValues);
