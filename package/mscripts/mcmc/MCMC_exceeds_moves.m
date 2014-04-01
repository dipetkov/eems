

function next = MCMC_exceeds_moves(mcmc,schedule)


next = 0;
type = mcmc.currType;


if (type<=5)
  if (schedule.paramtoupdate{type} == schedule.numparams{type})
    next = 1;
  end
else
  error('MCMC:Run','Unknown update type. Parameters not updated.')
end
