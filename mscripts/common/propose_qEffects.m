

function [proposal,pi1_pi0] = ...
    propose_qEffects(kernel,params,qVoronoi,mVoronoi,Sstruct,Mstruct,schedule)

%%%%%%%%%%
type = 2;%
%%%%%%%%%%

tile = schedule.paramtoupdate{type};

mEffcts = mVoronoi.mEffcts;
qEffcts = qVoronoi.qEffcts;
qEffcts(tile) = rnorm(qEffcts(tile),params.qEffctProposalS2);

proposal = struct('type',{type},'subtype',{1},'qEffcts',{qEffcts});

%% Constrain every tile effect to lie within range
if min( abs(qEffcts)<params.qEffctHalfInterval )
  mRates = realpow(10,mEffcts + params.mrateMu);
  qRates = realpow(10,qEffcts + params.qrateMu);
  [qValues,mValues] = ...
    average_rates(Mstruct,qRates,mRates,qVoronoi.qSeeds,mVoronoi.mSeeds,qVoronoi.Demes);
  proposal.kernel = resistance_kernel(Sstruct,Mstruct,mValues,qValues);
  pi1_pi0 = -sum(((qEffcts.^2)-(qVoronoi.qEffcts.^2))/(2*params.qrateS2));
else
  pi1_pi0 = -Inf;
end
