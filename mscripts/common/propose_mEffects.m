

function [proposal,pi1_pi0] = ...
    propose_mEffects(kernel,params,qVoronoi,mVoronoi,Sstruct,Mstruct,schedule)

%%%%%%%%%%
type = 6;%
%%%%%%%%%%

tile = schedule.paramtoupdate{type};

mEffcts = mVoronoi.mEffcts;
qEffcts = qVoronoi.qEffcts;
mEffcts(tile) = rnorm(mEffcts(tile),params.mEffctProposalS2);

proposal = struct('type',{type},'subtype',{1},'mEffcts',{mEffcts});
  
%% Constrain every tile effects to lie within range
if min( abs(mEffcts)<params.mEffctHalfInterval )
  mRates = realpow(10,mEffcts + params.mrateMu);
  qRates = realpow(10,qEffcts + params.qrateMu);
  [qValues,mValues] = ...
    average_rates(Mstruct,qRates,mRates,qVoronoi.qSeeds,mVoronoi.mSeeds,qVoronoi.Demes);
  proposal.kernel = resistance_kernel(Sstruct,Mstruct,mValues,qValues);
  pi1_pi0 = -sum(((mEffcts.^2)-(mVoronoi.mEffcts.^2))/(2*params.mrateS2));
else
  pi1_pi0 = -Inf;
end
