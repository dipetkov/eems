

function [mcmc,schedule] = index_schedule(Voronoi,params,mcmc,schedule)
%% It is necessary to update the schedule as the number of tiles changes %%


mtiles = Voronoi.mtiles;
schedule.numparams{2} = mtiles;
schedule.numparams{4} = mtiles;
