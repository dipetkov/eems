

function [Voronoi,params] = initial_values(Sstruct,Demes,habitat,opt)


p = Sstruct.nSites;
n = Sstruct.nIndiv;

if ~isfield(opt,'dfProposalS2')
  dfProposalS2 = sqrt(p);
else
  dfProposalS2 = opt.dfProposalS2;
end

%% Sample the initial number of tiles from zero-truncated negative
%% binomial distribution with high variance and mode away from zero
%% This ensures that independent runs will start with different 
%% Voronoi configurations (If always start with one or two tiles, 
%% then the initial configurations would be very similar)
mtiles = rnegbin(1,10,2/3);
%% Draw the Voronoi centers mSeeds uniformly within the habitat
mSeeds = runif_habitat(mtiles,habitat);
%% Assign migration rates to the Voronoi tiles
mrateS2 = rinvgam(1/2,1/2);
mrateMu = opt.mrateMuHalfInterval*(2*rand(1) - 1);
% On the log scale, the tile rates are mrateMu + mEffcts
mEffcts = init_effects(mtiles,opt.mEffctHalfInterval,mrateS2);
% The deviation of move proposals is scaled by the habitat range
mSeedsProposalS2 = opt.mSeedsProposalS2*[max(Demes(:,1))-min(Demes(:,1)),...
					 max(Demes(:,2))-min(Demes(:,2))];
Voronoi = struct('Demes',{Demes},...
		 'mtiles',{mtiles},...
                 'mSeeds',{mSeeds},...
                 'mEffcts',{mEffcts},...
                 'habitat',{habitat});
params = struct('mrateShape',{opt.mrateShape},...
		'mrateScale',{opt.mrateScale},...
                's2locShape',{opt.s2locShape},...
		's2locScale',{opt.s2locScale},...
                'negBiSize',{opt.negBiSize},...
		'negBiProb',{opt.negBiProb},...
		'mrateMu',{mrateMu},...
		'mrateS2',{mrateS2},...
		'mSeedsProposalS2',{mSeedsProposalS2},...
                'mEffctHalfInterval',{opt.mEffctHalfInterval},...
                'mrateMuHalfInterval',{opt.mrateMuHalfInterval},...
		'mEffctProposalS2',{opt.mEffctProposalS2},...
		'mrateMuProposalS2',{opt.mrateMuProposalS2},...
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
