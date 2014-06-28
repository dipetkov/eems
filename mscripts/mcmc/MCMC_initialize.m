

function [mcmc,schedule] = MCMC_initialize(Voronoi,params,opt)


numTypes = 5;
numThetas = params.numthetas;
okayMoves = cell(numTypes,1);
totalMoves = cell(numTypes,1);
dimTypes = ones(numTypes,1);
dimTypes(1) = numThetas;       %% The first type of update is a nuisance parameter  %%
dimTypes(end) = 2;             %% The final type of update is birth/death of a tile %%

for i=1:numTypes
  okayMoves{i} = zeros(dimTypes(i),1);
  totalMoves{i} = zeros(dimTypes(i),1);
end

mcmc = struct('currIter',{0},...
              'iterDone',{0},...
	      'iterStart',{0},...
	      'numBurn',{opt.numBurn},...
	      'numThin',{opt.numThin},...
	      'numIters',{opt.numIter},...
	      'isfinished',{0},...
	      'currType',{0},...
	      'numMoves',{0},...
	      'numTypes',{numTypes},...
	      'dimTypes',{dimTypes},...
	      'okayMoves',{okayMoves},...
              'totalMoves',{totalMoves});

mcmc.numSteps = floor((mcmc.numIters-mcmc.numBurn-1)/(mcmc.numThin+1));
mcmc.ii = mcmc.numBurn + 1 + (1:mcmc.numSteps)*(mcmc.numThin+1);
mcmc.kk = 1:mcmc.numSteps;

%% Initialize the tile and rate indices %%

numparams = cell(5,1);
numparams{1} = params.numthetas;
numparams{2} = Voronoi.ntiles;
numparams{3} = 1;
numparams{4} = Voronoi.ntiles;
numparams{5} = 1;

paramtoupdate = cell(5,1);
paramtoupdate{1} = 1;
paramtoupdate{2} = 1;
paramtoupdate{3} = 1;
paramtoupdate{4} = 1;
paramtoupdate{5} = 1;

schedule = struct('numparams',{numparams},...
                  'paramtoupdate',{paramtoupdate});
[mcmc,schedule] = index_schedule(Voronoi,params,mcmc,schedule);
