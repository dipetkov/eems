

function x = draw_rateMu(Min,Max,Mean,Var,s)
%% Proposal distribution for the mean log tile rate


if rand(1)<s
  %% long-range proposal %%
  x = rcauchy(Mean,1);
else
  %% local proposal %%
  x = rnorm(Mean,Var);
end
