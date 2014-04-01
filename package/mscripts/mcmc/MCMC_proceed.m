

function mcmc = MCMC_proceed(mcmc,proposal)


if (mcmc.currType~=proposal.type)
  error('method:MCMC','Failed: currType == proposal.type')
end

mcmc.currSubtype = proposal.subtype;
mcmc.totalMoves{mcmc.currType}(mcmc.currSubtype) = ...
    mcmc.totalMoves{mcmc.currType}(mcmc.currSubtype) + 1;
mcmc.iterStart = 0;
