

function mcmc = MCMC_fdisp(fid,mcmc)


fprintf(fid,'Acceptance proportions:\n');

for i=1:mcmc.numTypes
  a = mcmc.okayMoves{i};
  A = mcmc.totalMoves{i};
  fprintf(fid,'\t(%d/%d) = %.2f%% for type %d\n',...
          sum(a),sum(A),100*sum(a)/sum(A),i);
  if mcmc.dimTypes(i)>1
    for ij=1:mcmc.dimTypes(i)
      fprintf(fid,'\t\t(%d/%d) = %.2f%% for subtype %d\n',...
              a(ij),A(ij),100*a(ij)/A(ij),ij);
    end
  end
end
