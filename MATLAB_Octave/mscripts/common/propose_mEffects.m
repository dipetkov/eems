

function [proposal,pi1_pi0] = ...
    propose_mEffects(kernel,params,qVoronoi,mVoronoi,Data,Graph,schedule)

%%%%%%%%%%
type = 4;%
%%%%%%%%%%

tile = schedule.paramtoupdate{type};

mEffcts = mVoronoi.mEffcts;
qEffcts = qVoronoi.qEffcts;
mEffcts(tile) = rnorm(mEffcts(tile),params.mEffctProposalS2);
proposal = struct('type',{type},'subtype',{1},'mEffcts',{mEffcts});
  
%% Constrain every tile effects to lie within range
if min( abs(mEffcts)<params.mEffctHalfInterval )
  mRates = realpow(10,mEffcts + params.mrateMu);
  qRates = realpow(10,qEffcts); % qrateMu = 0.0;
  [qValues,mValues] = ...
    average_rates(Graph,qRates,mRates,qVoronoi.qSeeds,mVoronoi.mSeeds,qVoronoi.Demes);
  proposal.kernel = resistance_kernel(Data,Graph,mValues,qValues);
  pi1_pi0 = -sum(((mEffcts.^2)-(mVoronoi.mEffcts.^2))/(2*params.mrateS2));
else
  pi1_pi0 = -Inf;
end
