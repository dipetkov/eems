

%% Source directory where the *.m scripts are
sourcepath = '../';

%% Data directory where the input files are
%% Input filename; the input files are: 
%% dataset.coord: location coordinates, one individual per line
%% dataset.sites: matrix of alleles, one individual per line, two columns per microsatellite
%% dataset.dimns: 3x2 matrix tha specifies:
%%                      range of x coordinate
%%                      range of y coordinate
%%                      #samples #snps
%datapath = '~/package/data/sat-barrier-nIndiv150-nSites16-gridSize8x7';
datapath = '../data/sat-uniform-nIndiv150-nSites16-gridSize8x7';

%% Choose the size of the graph
%% The population graph is a regular triangular graph with size xPop-by-yPop
xPop = 12;
yPop = 8;

%% Output filename
simno = 1;  %% If you simulate several realizations of the Markov chain
mcmcpath = strcat(datapath,'-g',num2str(xPop),'x',num2str(yPop),'-simno',num2str(simno));

%% The input arguments sourcepath, datapath, mcmcpath, xPop, yPop have to be explicitly specified
%% There are optional input arguments (that have default values but can be tuned to improve convergence)
%% Here is a complete list:
%% numMCMCIter: number of MCMC iterations
%% numBurnIter: number of burn-in iterations
%% numThinIter: number of iterations to thin between two writing steps
%% mrateShape,mrateScale: hyperparameters for the cell effect variance omega_m^2
%% s2locShape,s2locShape: hyperparameters for the scale parameter sigma^2
%% negBiSize: number of failures for the Negative-Binomial prior on the number of Voronoi tiles
%% negBiProb: success probability for the Negative-Binomial prior on the number of Voronoi tiles
%% mEffctProposalS2:  proposal variance for the cell effects
%% mSeedsProposalS2:  proposal variance for the cell locations
%% mrateMuProposalS2: proposal variance for the overall migration log rate


addpath(sourcepath);
addpath(strcat(sourcepath,'/mscripts'));

MCMC_microsat(sourcepath,datapath,mcmcpath,xPop,yPop,...
	      'numMCMCIter',1200,'numBurnIter',600,'numThinIter',9);
