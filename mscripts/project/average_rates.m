

function [qValues,mValues] = average_rates(Graph,qRates,mRates,qSeeds,mSeeds,Demes)


%% The rate of an edge is the average of the rate of the two %%
%% tile it connets (possibly, it is one and the same tile)   %%
[temp,mColors] = pdist2(mSeeds,Demes,'euclidean');
mValues = 0.5*mRates(mColors(Graph.Vi)) + ...
          0.5*mRates(mColors(Graph.Vj));

[temp,qColors] = pdist2(qSeeds,Demes,'euclidean');
qValues = qRates(qColors);
