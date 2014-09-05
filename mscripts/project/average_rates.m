

function [qValues,mValues] = average_rates(Mstruct,qRates,mRates,qSeeds,mSeeds,Demes)


%% The rate of an edge is the average of the rate of the two %%
%% tile it connets (possibly, it is one and the same tile)   %%
[temp,mColors] = pdist2(mSeeds,Demes,'euclidean');
mValues = 0.5*mRates(mColors(Mstruct.Mi)) + ...
          0.5*mRates(mColors(Mstruct.Mj));

[temp,qColors] = pdist2(qSeeds,Demes,'euclidean');
qValues = qRates(qColors);
