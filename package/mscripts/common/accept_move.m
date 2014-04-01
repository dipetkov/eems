

function [Voronoi,kernel,params] = ...
    accept_move(Voronoi,kernel,params,proposal)


if (proposal.type==1)
  params = proposal.params;
elseif (proposal.type==2)
  Voronoi.Effcts = proposal.Effcts;
  kernel.X = proposal.X;
  kernel.XC = proposal.XC;
  kernel.ldDinv = proposal.ldDinv;
elseif (proposal.type==3)
  params.rateMu = proposal.rateMu;
  kernel.X = proposal.X;
  kernel.XC = proposal.XC;
  kernel.ldDinv = proposal.ldDinv;
elseif (proposal.type==4)
  Voronoi.Scoord = proposal.Scoord;
  kernel.X = proposal.X;
  kernel.XC = proposal.XC;
  kernel.ldDinv = proposal.ldDinv;
elseif (proposal.type==5)
  Voronoi.ntiles = proposal.ntiles;
  Voronoi.Scoord = proposal.Scoord;
  Voronoi.Effcts = proposal.Effcts;
  kernel.X = proposal.X;
  kernel.XC = proposal.XC;
  kernel.ldDinv = proposal.ldDinv;
end
