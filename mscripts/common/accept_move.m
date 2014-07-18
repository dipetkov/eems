

function [mVoronoi,kernel,params] = ...
    accept_move(mVoronoi,kernel,params,proposal)


if (proposal.type==1)
  params = proposal.params;
elseif (proposal.type==2)
  mVoronoi.mEffcts = proposal.mEffcts;
  kernel = proposal.kernel;
elseif (proposal.type==3)
  params.mrateMu = proposal.mrateMu;
  kernel = proposal.kernel;
elseif (proposal.type==4)
  mVoronoi.mSeeds = proposal.mSeeds;
  kernel = proposal.kernel;
elseif (proposal.type==5)
  mVoronoi.mtiles = proposal.mtiles;
  mVoronoi.mSeeds = proposal.mSeeds;
  mVoronoi.mEffcts = proposal.mEffcts;
  kernel = proposal.kernel;
end
