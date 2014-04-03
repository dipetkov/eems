

%% Source directory where the *.m scripts are
sourcepath = '~/package';

%% Data directory where the input files are
%% Input filename; the input files are: 
%% dataset.coord: location coordinates, one individual per line
%% dataset.sites: matrix of alleles, one individual per line, two columns per microsatellite
%% dataset.dimns: 3x2 matrix tha specifies:
%%                      range of x coordinate
%%                      range of y coordinate
%%                      #samples #snps
%datapath = '~/package/data/sat-tribars-nIndiv150-nSites16-gridSize8x7';
datapath = '~/package/data/sat-uniform-nIndiv150-nSites16-gridSize8x7';

%% Choose the size of the graph
%% The population graph is a regular triangular graph with size xPop-by-yPop
xPop = 12;
yPop = 8;

%% Output filename
simno = 1;  %% If you simulate several realizations of the Markov chain
mcmcpath = strcat(datapath,'-g',xPop,'x',yPop,'-simno',num2str(simno));


%% The input arguments sourcepath, datapath, mcmcpath, xPop, yPop have to be explicitly specified
%% There are optional input arguments (that have default values but can be tuned to improve convergence)
%% Here is a complete list:
%% numIter: number of MCMC iterations
%% numBurn: number of burn-in iterations
%% numThin: number of iterations to thin between two writing steps
%% effctS2: proposal variance for the cell effects
%% coordS2: proposal variance for the cell locations
%% log10S2: proposal variance for the overall migration log rate
%% ratesC: hyperparameter for the cell effect variance tau_m^2
%% s2locC: hyperparameter for the scale parameter sigma^2
%% negBiR: number of failures for the Negative-Binomial prior on the number of Voronoi tiles
%% negBiP: success probability for the Negative-Binomial prior on the number of Voronoi tiles


addpath(sourcepath);
addpath(strcat(sourcepath,'/mscripts'));

MCMC_microsats(sourcepath,datapath,mcmcpath,xPop,yPop,...
	       'numIter',1200,'numBurn',600,'numThin',9);
