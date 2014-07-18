

function [mcmc,schedule] = MCMC_change_update(mcmc,schedule)


type = mcmc.currType;

if (type<=5)
  schedule.paramtoupdate{type} = schedule.paramtoupdate{type} + 1;
  if (schedule.paramtoupdate{type} > schedule.numparams{type})
    schedule.paramtoupdate{type} = 1;
    mcmc.currType = type + 1;
  end
else
  error('MCMC:Run','Unknown update type. Parameters not updated.')
end

if (mcmc.currType>mcmc.numUpdateTypes)
  % End the iteration
  mcmc.currType = 1;
  mcmc.iterDone = 1;
end
