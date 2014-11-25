

function [mcmc,schedule] = index_schedule(qVoronoi,mVoronoi,params,mcmc,schedule)
%% It is necessary to update the schedule as the number of tiles changes %%


qtiles = qVoronoi.qtiles;
schedule.numparams{1} = qtiles;
schedule.numparams{2} = qtiles;
schedule.numparams{3} = qtiles;

mtiles = mVoronoi.mtiles;
schedule.numparams{4} = mtiles;
schedule.numparams{6} = mtiles;
schedule.numparams{7} = mtiles;
