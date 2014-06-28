

function [proposal,pi1_pi0] = ...
    propose_rateMu(kernel,params,Voronoi,Sstruct,Mstruct)

%%%%%%%%%%
type = 3;%
%%%%%%%%%%

Effcts = Voronoi.Effcts;
Scoord = Voronoi.Scoord;

rateMu = draw_rateMu(params.minRate,params.maxRate,...
                     params.rateMu,params.log10S2,...
    		     params.SmallWorld_s);

%% Constrain every effect to lie within range
if ( rateMu>params.minRate && rateMu<params.maxRate )
  mRates = realpow(10,Effcts + rateMu);
  mValues = average_rates(Mstruct,mRates,Scoord,Voronoi.Vcoord);
  proposal = resistance_kernel(Sstruct,Mstruct,mValues);
  proposal.type = type;
  proposal.subtype = 1;
  proposal.rateMu = rateMu;

  pi1_pi0 = 0;

else
  proposal = struct('type',{type},...
                    'subtype',{1},...
                    'rateMu',{rateMu});
  pi1_pi0 = -Inf;
end
