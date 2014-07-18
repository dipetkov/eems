

function [mcmc,schedule] = MCMC_initialize(Voronoi,params,opt)


numUpdateTypes = 5;
numThetas = params.numthetas;
okayMoves = cell(numUpdateTypes,1);
totalMoves = cell(numUpdateTypes,1);
dimUpdateTypes = ones(numUpdateTypes,1);
dimUpdateTypes(1) = numThetas;  %% The first type of update is a nuisance parameter  %%
dimUpdateTypes(end) = 2;        %% The final type of update is birth/death of a tile %%

for i=1:numUpdateTypes
  okayMoves{i} = zeros(dimUpdateTypes(i),1);
  totalMoves{i} = zeros(dimUpdateTypes(i),1);
end

mcmc = struct('currIter',{0},...
              'iterDone',{0},...
	      'iterStart',{0},...
	      'isfinished',{0},...
	      'currType',{0},...
	      'numMoves',{0},...
	      'numBurnIter',{opt.numBurnIter},...
	      'numThinIter',{opt.numThinIter},...
	      'numMCMCIter',{opt.numMCMCIter},...
	      'numUpdateTypes',{numUpdateTypes},...
	      'dimUpdateTypes',{dimUpdateTypes},...
	      'okayMoves',{okayMoves},...
              'totalMoves',{totalMoves});

mcmc.numSteps = floor((mcmc.numMCMCIter-mcmc.numBurnIter-1)/(mcmc.numThinIter+1));
mcmc.ii = mcmc.numBurnIter + 1 + (1:mcmc.numSteps)*(mcmc.numThinIter+1);
mcmc.kk = 1:mcmc.numSteps;

%% Initialize the tile and rate indices %%

numparams = cell(5,1);
numparams{1} = params.numthetas;
numparams{2} = Voronoi.mtiles;
numparams{3} = 1;
numparams{4} = Voronoi.mtiles;
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
