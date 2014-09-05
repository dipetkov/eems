

function Effects = init_effects(T,halfInterval,rateS2)
%% Prior distribution for the effects %%


Effects = rtnorm(T,-halfInterval,halfInterval,0,rateS2);
