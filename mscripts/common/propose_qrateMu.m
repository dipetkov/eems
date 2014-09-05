

function [proposal,pi1_pi0] = ...
    propose_qrateMu(kernel,params,qVoronoi,mVoronoi,Sstruct,Mstruct)

%%%%%%%%%%
type = 3;%
%%%%%%%%%%

mrateMu = params.mrateMu;
qrateMu = rnorm(params.qrateMu,params.qrateMuProposalS2);

proposal = struct('type',{type},'subtype',{1},'qrateMu',{qrateMu});

%% Constrain every tile effect to lie within range
if ( abs(qrateMu)<params.qrateMuHalfInterval )
  mRates = realpow(10,mVoronoi.mEffcts + mrateMu);
  qRates = realpow(10,qVoronoi.qEffcts + qrateMu);
  [qValues,mValues] = ...
    average_rates(Mstruct,qRates,mRates,qVoronoi.qSeeds,mVoronoi.mSeeds,qVoronoi.Demes);
  proposal.kernel = resistance_kernel(Sstruct,Mstruct,mValues,qValues);
  pi1_pi0 = 0;
else
  pi1_pi0 = -Inf;
end
