

function [proposal,pi1_pi0] = ...
    move_qVoronoi(kernel,params,qVoronoi,mVoronoi,Sstruct,Mstruct,schedule)

%%%%%%%%%%
type = 4;%
%%%%%%%%%%

tile = schedule.paramtoupdate{type};

mSeeds = mVoronoi.mSeeds;
qSeeds = qVoronoi.qSeeds;
qSeeds(tile,:) = rnorm(qSeeds(tile,:),params.qSeedsProposalS2);

proposal = struct('type',{type},'subtype',{1},'qSeeds',{qSeeds});

%% Contain every center to lie within the habitat
if min(is_in_habitat(qVoronoi.habitat,qSeeds))
  mRates = realpow(10,mVoronoi.mEffcts + params.mrateMu);
  qRates = realpow(10,qVoronoi.qEffcts + params.qrateMu);
  [qValues,mValues] = ...
    average_rates(Mstruct,qRates,mRates,qSeeds,mSeeds,qVoronoi.Demes);
  proposal.kernel = resistance_kernel(Sstruct,Mstruct,mValues,qValues);
  pi1_pi0 = 0;
else
  pi1_pi0 = -Inf;
end
