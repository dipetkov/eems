

function [qVoronoi,mVoronoi,kernel,params] = ...
    accept_move(qVoronoi,mVoronoi,kernel,params,proposal)


if (proposal.type==1)
  qVoronoi.qEffcts = proposal.qEffcts;
  kernel = proposal.kernel;
elseif (proposal.type==2)
  qVoronoi.qSeeds = proposal.qSeeds;
  kernel = proposal.kernel;
elseif (proposal.type==3)
  qVoronoi.qtiles = proposal.qtiles;
  qVoronoi.qSeeds = proposal.qSeeds;
  qVoronoi.qEffcts = proposal.qEffcts;
  kernel = proposal.kernel;
elseif (proposal.type==4)
  mVoronoi.mEffcts = proposal.mEffcts;
  kernel = proposal.kernel;
elseif (proposal.type==5)
  params.mrateMu = proposal.mrateMu;
  kernel = proposal.kernel;
elseif (proposal.type==6)
  mVoronoi.mSeeds = proposal.mSeeds;
  kernel = proposal.kernel;
elseif (proposal.type==7)
  mVoronoi.mtiles = proposal.mtiles;
  mVoronoi.mSeeds = proposal.mSeeds;
  mVoronoi.mEffcts = proposal.mEffcts;
  kernel = proposal.kernel;
else
  params = proposal.params;
end
