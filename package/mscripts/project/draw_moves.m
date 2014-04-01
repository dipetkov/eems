

function x = draw_moves(habitat,Mean,Var,s)
%% Proposal distribution for tile moves %%


if runif(1)<s
  %% long-range proposal %%
  x = rcauchy(Mean,1);
else
  %% local proposal %%
  x = rnorm(Mean,Var);
end
