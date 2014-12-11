

sourcepath = '../';
addpath(strcat(sourcepath,'/mscripts'));
datapath = '/Users/desislava/testeems/data/barrier-schemeX-nIndiv310-s12x8-u4Nm1-L3000';
xDemes = 13;
yDemes = 7;
simno = 0;
mcmcpath = strcat(datapath,'-g',num2str(xDemes),'x',num2str(yDemes),'-simno',num2str(simno));


MCMC_haploid(sourcepath,datapath,mcmcpath,xDemes,yDemes,...
	     'numMCMCIter',12000,'numBurnIter',6000,'numThinIter',99);
