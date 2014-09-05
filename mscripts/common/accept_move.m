

function [qVoronoi,mVoronoi,kernel,params] = ...
    accept_move(qVoronoi,mVoronoi,kernel,params,proposal)


if (proposal.type==1)
  params = proposal.params;
elseif (proposal.type==2)
  qVoronoi.qEffcts = proposal.qEffcts;
  kernel = proposal.kernel;
elseif (proposal.type==3)
  params.qrateMu = proposal.qrateMu;
  kernel = proposal.kernel;
elseif (proposal.type==4)
  qVoronoi.qSeeds = proposal.qSeeds;
  kernel = proposal.kernel;
elseif (proposal.type==5)
  qVoronoi.qtiles = proposal.qtiles;
  qVoronoi.qSeeds = proposal.qSeeds;
  qVoronoi.qEffcts = proposal.qEffcts;
  kernel = proposal.kernel;
elseif (proposal.type==6)
  mVoronoi.mEffcts = proposal.mEffcts;
  kernel = proposal.kernel;
elseif (proposal.type==7)
  params.mrateMu = proposal.mrateMu;
  kernel = proposal.kernel;
elseif (proposal.type==8)
  mVoronoi.mSeeds = proposal.mSeeds;
  kernel = proposal.kernel;
elseif (proposal.type==9)
  mVoronoi.mtiles = proposal.mtiles;
  mVoronoi.mSeeds = proposal.mSeeds;
  mVoronoi.mEffcts = proposal.mEffcts;
  kernel = proposal.kernel;
end
