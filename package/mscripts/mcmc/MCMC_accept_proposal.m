

function mcmc = MCMC_accept_proposal(mcmc)


mcmc.okayMoves{mcmc.currType}(mcmc.currSubtype) = ...
    mcmc.okayMoves{mcmc.currType}(mcmc.currSubtype) + 1;

