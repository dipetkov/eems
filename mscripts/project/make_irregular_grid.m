

function [Mstruct,Demes,Edges] = make_irregular_grid(datapath,inDemes,Mij,Demes,Edges)
%% Remove vertices and edges that fall outside an irregularly shaped habitat


nDemes = length(inDemes);
nEdges = nrow(Mij);
inEdges = ones(nEdges,1);

for alpha = 1:nDemes
  if ~inDemes(alpha)
    inPair = sum(Mij==alpha,2);
    inEdges(inPair>0) = 0;
  end
end

%% This is the number of demes to keep
%% It is important however to make sure
%% that the resulting graph is connected
nEdges = sum(inEdges);
nDemes = sum(inDemes);

if nDemes~=nrow(inDemes)
  %% It is necessary to re-index the demes
  C = 1:nDemes;
  indxDemes = inDemes;
  indxDemes(inDemes==1) = C;
  %% Remove all rows of Mij that map to an outside deme
  Mij = Mij(inEdges==1,:);
  Mij = indxDemes(Mij);
  %% Re-index both the vertices and the edges
  %% (The edges are only needed for plotting)
  Demes = Demes(inDemes==1,:);
  Edges = Edges(inDemes==1,:);
  allDemes = length(inDemes);
  indxDemes = [indxDemes;0];
  Edges(Edges==0) = allDemes+1;
  Edges = indxDemes(Edges);
end

%% Check that the population graph is connected
M = sparse(Mij(:,1),Mij(:,2),1) + speye(nDemes);
[p,q,r,s] = dmperm(M);
if length(r)>2
  error('The population graph is not connected.')
end

Mstruct = struct('Mi',{Mij(:,1)},...
                 'Mj',{Mij(:,2)},...
                 'nDemes',{nDemes},...
		 'nEdges',{nEdges});
