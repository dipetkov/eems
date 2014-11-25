

function [mcmc,schedule] = MCMC_change_update(mcmc,schedule)


type = mcmc.currType;
schedule.paramtoupdate{type} = schedule.paramtoupdate{type} + 1;
if (schedule.paramtoupdate{type} > schedule.numparams{type})
  schedule.paramtoupdate{type} = 1;
  mcmc.currType = type + 1;
end
if (mcmc.currType>mcmc.numUpdateTypes)
  mcmc.currType = 1;
  mcmc.iterDone = 1;
end
