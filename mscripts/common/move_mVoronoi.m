

function [proposal,pi1_pi0] = ...
    move_mVoronoi(kernel,params,mVoronoi,Sstruct,Mstruct,schedule)

%%%%%%%%%%
type = 4;%
%%%%%%%%%%

tile = schedule.paramtoupdate{type};

mSeeds = mVoronoi.mSeeds;
mSeeds(tile,:) = rnorm(mSeeds(tile,:),params.mSeedsProposalS2);

proposal = struct('type',{type},'subtype',{1},'mSeeds',{mSeeds});

% Contain every center to lie within the habitat
if min(is_in_habitat(mVoronoi.habitat,mSeeds))
  mRates = realpow(10,params.mrateMu + mVoronoi.mEffcts);
  mValues = average_rates(Mstruct,mRates,mSeeds,mVoronoi.Demes);
  proposal.kernel = resistance_kernel(Sstruct,Mstruct,mValues);
  pi1_pi0 = 0;
else
  pi1_pi0 = -Inf;
end
