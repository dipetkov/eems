

function mValues = average_rates(Mstruct,mRates,Scoord,Vcoord)
%% Scoord = Voronoi centers %%
%% Vcoord = graph vertices  %%


%% The rate of an edge is the average of the rate of the two %%
%% tile it connets (possibly, it is one and the same tile)   %%
[temp,Colors] = pdist2(Scoord,Vcoord,'euclidean');
mValues = 0.5*mRates(Colors(Mstruct.Mi)) + ...
          0.5*mRates(Colors(Mstruct.Mj));
