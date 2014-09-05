

function [proposal,pi1_pi0] = ...
    propose_mrateMu(kernel,params,qVoronoi,mVoronoi,Sstruct,Mstruct)

%%%%%%%%%%
type = 7;%
%%%%%%%%%%

qrateMu = params.qrateMu;
mrateMu = rnorm(params.mrateMu,params.mrateMuProposalS2);

proposal = struct('type',{type},'subtype',{1},'mrateMu',{mrateMu});

%% Constrain every tile effect to lie within range
if ( abs(mrateMu)<params.mrateMuHalfInterval )
  mRates = realpow(10,mVoronoi.mEffcts + mrateMu);
  qRates = realpow(10,qVoronoi.qEffcts + qrateMu);
  [qValues,mValues] = ...
    average_rates(Mstruct,qRates,mRates,qVoronoi.qSeeds,mVoronoi.mSeeds,mVoronoi.Demes);
  proposal.kernel = resistance_kernel(Sstruct,Mstruct,mValues,qValues);
  pi1_pi0 = 0;
else
  pi1_pi0 = -Inf;
end
