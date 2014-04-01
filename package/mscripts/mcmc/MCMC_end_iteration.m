

function mcmc = MCMC_end_iteration(mcmc,energy)


if ~mod(mcmc.currIter,1000)
  fprintf(2,'Ending iteration %d with log likelihood = %7.5f\n',mcmc.currIter,energy);
  fprintf(2,'   and acceptance proportions\n');
  for i=1:mcmc.numTypes
    a = mcmc.okayMoves{i};
    A = mcmc.totalMoves{i};
    fprintf(2,'\t(%d/%d) = %.2f%% for type %d\n',...
            sum(a),sum(A),100*sum(a)/sum(A),i);
    if mcmc.dimTypes(i)>1
      for ij=1:mcmc.dimTypes(i)
        fprintf(2,'\t\t(%d/%d) = %.2f%% for subtype %d\n',...
                a(ij),A(ij),100*a(ij)/A(ij),ij);
      end
    end
  end
end

if (mcmc.currIter==mcmc.numIters)
  mcmc.isfinished = 1;
end
