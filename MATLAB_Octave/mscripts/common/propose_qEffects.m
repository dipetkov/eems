

function [proposal,pi1_pi0] = ...
    propose_qEffects(kernel,params,qVoronoi,mVoronoi,Data,Graph,schedule)

%%%%%%%%%%
type = 1;%
%%%%%%%%%%

tile = schedule.paramtoupdate{type};

mEffcts = mVoronoi.mEffcts;
qEffcts = qVoronoi.qEffcts;
qEffcts(tile) = rnorm(qEffcts(tile),params.qEffctProposalS2);
proposal = struct('type',{type},'subtype',{1},'qEffcts',{qEffcts});

%% Constrain every tile effect to lie within range
if min( abs(qEffcts)<params.qEffctHalfInterval )
  mRates = realpow(10,mEffcts + params.mrateMu);
  qRates = realpow(10,qEffcts); % qrateMu = 0.0;
  [qValues,mValues] = ...
    average_rates(Graph,qRates,mRates,qVoronoi.qSeeds,mVoronoi.mSeeds,qVoronoi.Demes);
  proposal.kernel = resistance_kernel(Data,Graph,mValues,qValues);
  pi1_pi0 = -sum(((qEffcts.^2)-(qVoronoi.qEffcts.^2))/(2*params.qrateS2));
else
  pi1_pi0 = -Inf;
end
