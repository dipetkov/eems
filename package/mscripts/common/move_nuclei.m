

function [proposal,pi1_pi0] = ...
    move_nuclei(kernel,params,Voronoi,Sstruct,Mstruct,schedule)

%%%%%%%%%%
type = 4;%
%%%%%%%%%%

Effcts = Voronoi.Effcts;
Scoord = Voronoi.Scoord;
rateMu = params.rateMu;

tile = schedule.paramtoupdate{type};
Scoord(tile,:) = draw_moves(Voronoi.habitat,Scoord(tile,:),...
			    params.coordS2,params.SmallWorld_s);

% Contain every center to lie within the habitat
if min(is_in_habitat(Voronoi.habitat,Scoord))
  mRates = realpow(10,params.rateMu + Effcts);
  mValues = average_rates(Mstruct,mRates,Scoord,Voronoi.Vcoord);
  proposal = resistance_kernel(Sstruct,Mstruct,mValues);
  proposal.type = type;
  proposal.subtype = 1;
  proposal.Scoord = Scoord;

  pi1_pi0 = 0;

else
  proposal = struct('type',{type},'subtype',{1},...
                    'Scoord',{Scoord});
  pi1_pi0 = -Inf;
end
