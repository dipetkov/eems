

function mValues = average_rates(Mstruct,mRates,mSeeds,Demes)
%% mSeeds = Voronoi centers %%
%% Demes = graph vertices  %%


%% The rate of an edge is the average of the rate of the two %%
%% tile it connets (possibly, it is one and the same tile)   %%
[temp,Colors] = pdist2(mSeeds,Demes,'euclidean');
mValues = 0.5*mRates(Colors(Mstruct.Mi)) + ...
          0.5*mRates(Colors(Mstruct.Mj));
