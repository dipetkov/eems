

function [mcmc,schedule] = MCMC_initialize(qVoronoi,mVoronoi,params,opt)


numUpdateTypes = 9;
numThetas = params.numthetas;
okayMoves = cell(numUpdateTypes,1);
totalMoves = cell(numUpdateTypes,1);
dimUpdateTypes = ones(numUpdateTypes,1);
dimUpdateTypes(1) = numThetas;     %% This update is a nuisance parameter   %%
dimUpdateTypes(5) = 2;             %% This update is birth/death of a qtile %%
dimUpdateTypes(9) = 2;             %% This update is birth/death of a mtile %%

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

qtiles = qVoronoi.qtiles;
mtiles = mVoronoi.mtiles;

numparams = cell(9,1);
numparams{1} = params.numthetas;
numparams{2} = qtiles;
numparams{3} = 1;
numparams{4} = qtiles;
numparams{5} = qtiles;

numparams{6} = mtiles;
numparams{7} = 1;
numparams{8} = mtiles;
numparams{9} = mtiles;

paramtoupdate = cell(9,1);
paramtoupdate{1} = 1;
paramtoupdate{2} = 1;
paramtoupdate{3} = 1;
paramtoupdate{4} = 1;
paramtoupdate{5} = 1;
paramtoupdate{6} = 1;
paramtoupdate{7} = 1;
paramtoupdate{8} = 1;
paramtoupdate{9} = 1;

schedule = struct('numparams',{numparams},...
                  'paramtoupdate',{paramtoupdate});
[mcmc,schedule] = index_schedule(qVoronoi,mVoronoi,params,mcmc,schedule);
