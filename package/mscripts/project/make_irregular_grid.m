

function [Mstruct,Vcoord,Vedges] = ...
    make_irregular_grid(datapath,inPops,Mij,Vcoord,Vedges)
%% Remove vertices and edges that fall outside an irregularly %%
%% shaped habitat                                             %%


nPop = length(inPops);
nEdges = nrow(Mij);
inEdges = ones(nEdges,1);

for alpha = 1:nPop
  if ~inPops(alpha)
    inPair = sum(Mij==alpha,2);
    inEdges(inPair>0) = 0;
  end
end

%% This is the number of demes to keep
%% It is important however to make sure
%% that the resulting graph is connected
nEdges = sum(inEdges);
nPop = sum(inPops);

if nPop~=nrow(inPops)
  %% It is necessary to re-index the demes
  C = 1:nPop;
  indxPops = inPops;
  indxPops(inPops==1) = C;
  %% Remove all rows of Mij that map to an outside deme
  Mij = Mij(inEdges==1,:);
  Mij = indxPops(Mij);
  %% Re-index both the vertices and the edges
  %% (The edges are only needed for plotting)
  Vcoord = Vcoord(inPops==1,:);
  Vedges = Vedges(inPops==1,:);
  allPop = length(inPops);
  indxPops = [indxPops;0];
  Vedges(Vedges==0) = allPop+1;
  Vedges = indxPops(Vedges);
end

%% Check that the population graph is connected
M = sparse(Mij(:,1),Mij(:,2),1) + speye(nPop);
[p,q,r,s] = dmperm(M);
if length(r)>2
  error('The population graph is not connected.')
end

Mstruct = struct('Mi',{Mij(:,1)},...
                 'Mj',{Mij(:,2)},...
                 'nDemes',{nPop},...
		 'nEdges',{nEdges});
