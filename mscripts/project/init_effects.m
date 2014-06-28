

function Effcts = init_effects(T,absDiff,rateS2)
%% Prior distribution for the cell effects %%


lob = -absDiff;
upb =  absDiff;
Effcts = rtnorm(T,lob,upb,0,rateS2);
