

function [proposal,pi1_pi0] = ...
    propose_mrateMu(kernel,params,mVoronoi,Sstruct,Mstruct)

%%%%%%%%%%
type = 3;%
%%%%%%%%%%

mrateMu = rnorm(params.mrateMu,params.mrateMuProposalS2);

proposal = struct('type',{type},'subtype',{1},'mrateMu',{mrateMu});

%% Constrain every effect to lie within range
if ( abs(mrateMu)<params.mrateMuHalfInterval )
  mRates = realpow(10,mVoronoi.mEffcts + mrateMu);
  mValues = average_rates(Mstruct,mRates,mVoronoi.mSeeds,mVoronoi.Demes);
  proposal.kernel = resistance_kernel(Sstruct,Mstruct,mValues);
  pi1_pi0 = 0;
else
  pi1_pi0 = -Inf;
end
