

function [proposal,pi1_pi0] = ...
    propose_mEffects(kernel,params,mVoronoi,Sstruct,Mstruct,schedule)

%%%%%%%%%%
type = 2;%
%%%%%%%%%%

tile = schedule.paramtoupdate{type};

mEffcts = mVoronoi.mEffcts;
mEffcts(tile) = rnorm(mEffcts(tile),params.mEffctProposalS2);

proposal = struct('type',{type},'subtype',{1},'mEffcts',{mEffcts});
  
%% Constrain every effect to lie within range
if min( abs(mEffcts)<params.mEffctHalfInterval)
  mRates = realpow(10,params.mrateMu + mEffcts);
  mValues = average_rates(Mstruct,mRates,mVoronoi.mSeeds,mVoronoi.Demes);
  proposal.kernel = resistance_kernel(Sstruct,Mstruct,mValues);
  pi1_pi0 = -sum((mEffcts.^2-mVoronoi.mEffcts.^2)/(2*params.mrateS2));
else
  pi1_pi0 = -Inf;
end
