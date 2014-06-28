

function [mcmc,schedule] = index_schedule(Voronoi,params,mcmc,schedule)
%% It is necessary to update the schedule as the number of tiles changes %%


T = Voronoi.ntiles;
schedule.numparams{2} = T;
schedule.numparams{4} = T;
