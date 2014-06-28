

function [proposal,pi1_pi0] = ...
    propose_effects(kernel,params,Voronoi,Sstruct,Mstruct,schedule)

%%%%%%%%%%
type = 2;%
%%%%%%%%%%

Effcts = Voronoi.Effcts;
Scoord = Voronoi.Scoord;
rateMu = params.rateMu;

tile = schedule.paramtoupdate{type};
Effcts(tile) = draw_effects(params.absDiff,Effcts(tile),...
			    params.effctS2,params.SmallWorld_s);
  
%% Constrain every effect to lie within range
if min( abs(Effcts)<params.absDiff )
  mRates = realpow(10,rateMu + Effcts);
  mValues = average_rates(Mstruct,mRates,Scoord,Voronoi.Vcoord);
  proposal = resistance_kernel(Sstruct,Mstruct,mValues);
  proposal.type = type;
  proposal.subtype = 1;
  proposal.Effcts = Effcts;

  pi1_pi0 = -(proposal.Effcts(tile)^2-Voronoi.Effcts(tile)^2)/(2*params.rateS2);

else
  proposal = struct('type',{type},...
                    'subtype',{1},...
                    'Effcts',{Effcts});
  pi1_pi0 = -Inf;
end
