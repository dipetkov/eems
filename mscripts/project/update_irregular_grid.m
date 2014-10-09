

function [Graph,Demes,Edges] = ...
	 update_irregular_grid(datapath,inDeme,Demes,Edges,Pairs)
%% Remove vertices and edges that fall outside an irregularly %%
%% shaped habitat                                             %%


nDemes = nrow(Demes);
nEdges = nrow(Pairs);

if nDemes ~= sum(inDeme)
  Ci = 1:sum(inDeme);
  newIndex = inDeme;
  newIndex(inDeme==1) = Ci;
  %% It is necessary to reindex the demes
  [Demes,Edges,Pairs] = reindex_demes(inDeme,newIndex,Demes,Edges,Pairs);
  nDemes = nrow(Demes);
  nEdges = nrow(Pairs);
end

%% Check that the population graph is connected
G = sparse(Pairs(:,1),Pairs(:,2),1) + speye(nDemes);
[p,q,r,s] = dmperm(G);
if length(r)>2
  error('The population graph is not connected.')
end

Graph = struct('Vi',{Pairs(:,1)},...
               'Vj',{Pairs(:,2)},...
               'nDemes',{nDemes},...
	       'nEdges',{nEdges});
