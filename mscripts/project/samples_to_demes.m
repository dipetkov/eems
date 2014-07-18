

function [Juniq,Jinvpt,Jindex] = samples_to_demes(Samples,Demes)
%% Demes: the locations of the demes
%% Samples: the sampling locations 


%% Juniq: indices of sampled demes (those with at least one observation)
%% Jinvpt: map from demes to samples
%% Index: map from samples to demes
%% so that
%% Juniq(Jinvpt)==Jindex
[tempi,Jindex] = pdist2(Demes,Samples,'euclidean');
[Juniq,tempi,Jinvpt] = unique(Jindex);
