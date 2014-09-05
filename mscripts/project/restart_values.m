

function [qVoronoi,mVoronoi,params] = restart_values(Sstruct,Demes,habitat,opt,mcmcpath)


p = Sstruct.nSites;
n = Sstruct.nIndiv;

if ~isfield(opt,'dfS2')
  dfProposalS2 = sqrt(p);
else
  dfProposalS2 = opt.dfS2;
end

mcmcmtiles = dlmread(strcat(mcmcpath,'.mcmcmtiles'));
mcmcmhyper = dlmread(strcat(mcmcpath,'.mcmcmhyper'));

nIters = length(mcmcmtiles);
mtiles = mcmcmtiles(nIters);
mhyper = mcmcmhyper(nIters,:);

mcmcxCoord = fopen(strcat(mcmcpath,'.mcmcxcoord'));
mcmcyCoord = fopen(strcat(mcmcpath,'.mcmcycoord'));
mcmcmRates = fopen(strcat(mcmcpath,'.mcmcmrates'));

for line=1:nIters-1
  tempi = fgetl(mcmcxCoord);
  tempi = fgetl(mcmcyCoord);
  tempi = fgetl(mcmcmRates);
end

xCoord = str2num(fgetl(mcmcxCoord))';
yCoord = str2num(fgetl(mcmcyCoord))';
mRates = str2num(fgetl(mcmcmRates))';

fclose(mcmcxCoord);
fclose(mcmcyCoord);
fclose(mcmcmRates);

%% Sample the initial number of tiles from zero-truncated negative
%% binomial distribution with high variance and mode away from zero
%% This ensures that independent runs will start with different 
%% Voronoi configurations (If always start with one or two tiles, 
%% then the initial configurations would be very similar)
qtiles = rnegbin(1,10,2/3);
%% Draw the Voronoi centers Coord uniformly within the habitat
qSeeds = runif_habitat(qtiles,habitat);
mSeeds = [xCoord,yCoord];
%% Assign migration rates to the Voronoi tiles
mrateS2 = mhyper(2);
qrateS2 = rinvgam(1/2,1/2);
mrateMu = mhyper(1);
qrateMu = opt.qrateMuHalfInterval*(2*rand(1) - 1);
mEffcts = log10(mRates) - mrateMu;
qEffcts = init_effects(qtiles,opt.qEffctHalfInterval,qrateS2);
% The deviation of move proposals is scaled by the habitat range
mSeedsProposalS2 = opt.mSeedsProposalS2*[max(Demes(:,1))-min(Demes(:,1)),...
					 max(Demes(:,2))-min(Demes(:,2))];
qSeedsProposalS2 = opt.qSeedsProposalS2*[max(Demes(:,1))-min(Demes(:,1)),...
					 max(Demes(:,2))-min(Demes(:,2))];
qVoronoi = struct('Demes',{Demes},...
		  'qtiles',{qtiles},...
                  'qSeeds',{qSeeds},...
                  'qEffcts',{qEffcts},...
                  'habitat',{habitat});
mVoronoi = struct('Demes',{Demes},...
		  'mtiles',{mtiles},...
                  'mSeeds',{mSeeds},...
                  'mEffcts',{mEffcts},...
                  'habitat',{habitat});
params = struct('mrateShape',{opt.mrateShape},...
		'mrateScale',{opt.mrateScale},...
                'qrateShape',{opt.qrateShape},...
		'qrateScale',{opt.qrateScale},...
                's2locShape',{opt.s2locShape},...
		's2locScale',{opt.s2locScale},...
                'negBiSize',{opt.negBiSize},...
		'negBiProb',{opt.negBiProb},...
		'mrateMu',{mrateMu},...
		'mrateS2',{mrateS2},...
		'qrateMu',{qrateMu},...
		'qrateS2',{qrateS2},...
		'mEffctHalfInterval',{opt.mEffctHalfInterval},...
		'qEffctHalfInterval',{opt.qEffctHalfInterval},...
		'mrateMuHalfInterval',{opt.mrateMuHalfInterval},...
		'qrateMuHalfInterval',{opt.qrateMuHalfInterval},...
		'mrateMuProposalS2',{opt.mrateMuProposalS2},...
		'qrateMuProposalS2',{opt.qrateMuProposalS2},...
		'mEffctProposalS2',{opt.mEffctProposalS2},...
		'qEffctProposalS2',{opt.qEffctProposalS2},...
		'mSeedsProposalS2',{mSeedsProposalS2},...
		'qSeedsProposalS2',{qSeedsProposalS2},...
		'dfProposalS2',{dfProposalS2});

if Sstruct.microsat
  params.df = p;
  params.dflob = p - 1;
  params.dfmin = p - 1;
  params.dfmax = p + 1;
  params.dfupb = p + 1;
else
  params.df = n;
  params.dflob = n - 1;
  params.dfmin = n - 1;
  params.dfmax = n + 1;
  params.dfupb = p;
end
