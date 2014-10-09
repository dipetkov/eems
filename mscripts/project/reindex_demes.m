

function [Demes,Edges,Pairs] = reindex_demes(inDeme,newIndex,Demes,Edges,Pairs)


nDemes = nrow(Demes);
inDeme = (inDeme>0);
inEdge = inDeme(Pairs(:,1)).*inDeme(Pairs(:,2));
nDemes2 = sum(inDeme);
nEdges2 = sum(inEdge);

%% Remove all rows of Pairs that map to an outside deme
Pairs = Pairs(inEdge==1,:);
Pairs = newIndex(Pairs);
%% Re-index both the vertices and the edges
%% (The edges are only needed for plotting)
Demes = Demes(inDeme==1,:);
Edges = Edges(inDeme==1,:);
newIndex = [newIndex; 0];
Edges(Edges==0)=nDemes+1;
Edges = newIndex(Edges);
