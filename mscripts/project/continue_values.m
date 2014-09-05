

function [qVoronoi,mVoronoi,params] = continue_values(Sstruct,Demes,habitat,opt,mcmcpath)


p = Sstruct.nSites;
n = Sstruct.nIndiv;

if ~isfield(opt,'dfS2')
  dfProposalS2 = sqrt(p);
else
  dfProposalS2 = opt.dfS2;
end

mcmcthetas = dlmread(strcat(mcmcpath,'.mcmcthetas'));
mcmcmtiles = dlmread(strcat(mcmcpath,'.mcmcmtiles'));
mcmcqtiles = dlmread(strcat(mcmcpath,'.mcmcqtiles'));
mcmcmhyper = dlmread(strcat(mcmcpath,'.mcmcmhyper'));
mcmcqhyper = dlmread(strcat(mcmcpath,'.mcmcqhyper'));

nIters = length(mcmcmtiles);
mtiles = mcmcmtiles(nIters);
qtiles = mcmcqtiles(nIters);
mhyper = mcmcmhyper(nIters,:);
qhyper = mcmcqhyper(nIters,:);
thetas = mcmcthetas(nIters,:);

mcmcxCoord = fopen(strcat(mcmcpath,'.mcmcxcoord'));
mcmcyCoord = fopen(strcat(mcmcpath,'.mcmcycoord'));
mcmcwCoord = fopen(strcat(mcmcpath,'.mcmcwcoord'));
mcmczCoord = fopen(strcat(mcmcpath,'.mcmczcoord'));
mcmcmRates = fopen(strcat(mcmcpath,'.mcmcmrates'));
mcmcqRates = fopen(strcat(mcmcpath,'.mcmcqrates'));

for line=1:nIters-1
  tempi = fgetl(mcmcxCoord);
  tempi = fgetl(mcmcyCoord);
  tempi = fgetl(mcmcmRates);
  tempi = fgetl(mcmcwCoord);
  tempi = fgetl(mcmczCoord);
  tempi = fgetl(mcmcqRates);
end

xCoord = str2num(fgetl(mcmcxCoord))';
yCoord = str2num(fgetl(mcmcyCoord))';
wCoord = str2num(fgetl(mcmcwCoord))';
zCoord = str2num(fgetl(mcmczCoord))';
mRates = str2num(fgetl(mcmcmRates))';
qRates = str2num(fgetl(mcmcqRates))';

fclose(mcmcxCoord);
fclose(mcmcyCoord);
fclose(mcmcwCoord);
fclose(mcmczCoord);
fclose(mcmcmRates);
fclose(mcmcqRates);

%% Sample the initial number of tiles from zero-truncated negative
%% binomial distribution with high variance and mode away from zero
%% This ensures that independent runs will start with different 
%% Voronoi configurations (If always start with one or two tiles, 
%% then the initial configurations would be very similar)
%% Draw the Voronoi centers Coord uniformly within the habitat
qSeeds = [wCoord,zCoord];
mSeeds = [xCoord,yCoord];
%% Assign migration rates to the Voronoi tiles
mrateS2 = mhyper(2);
qrateS2 = qhyper(2);
mrateMu = mhyper(1);
qrateMu = qhyper(1);
mEffcts = log10(mRates) - mrateMu;
qEffcts = log10(qRates) - qrateMu;
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
  params.df = thetas(2);
  params.dflob = n - 1;
  params.dfmin = n - 1;
  params.dfmax = p;
  params.dfupb = p;
end
