

function mcmc = MCMC_start_iteration(mcmc)


if (mcmc.currIter==0)
  mcmc.currType = 1;
end

mcmc.currIter = mcmc.currIter+1;
mcmc.iterDone = 0;

if mod(mcmc.currIter,1000)==1
  fprintf(2,'Starting iteration %d of %d\n',mcmc.currIter,mcmc.numMCMCIter);
end
