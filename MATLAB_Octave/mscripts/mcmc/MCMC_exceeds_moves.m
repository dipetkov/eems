

function next = MCMC_exceeds_moves(mcmc,schedule)


next = 0;
type = mcmc.currType;

if (schedule.paramtoupdate{type} == schedule.numparams{type})
  next = 1;
end
