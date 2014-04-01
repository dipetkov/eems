

function x = draw_effects(Abs,Mean,Var,s)
%% Proposal distribution for the cell effects %%


if runif(1)<s
  %% long-range proposal %%
  x = rcauchy(Mean,1);
else
  %% local proposal %%
  x = rnorm(Mean,Var);
end
