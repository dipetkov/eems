

function [Graph,Demes,Edges,Pairs,habitat] = make_triangular_grid(datapath,xDemes,yDemes)
%% Construct three data structures that describe the triangular grid       %%
%% Demes: list of vertices (one vertex per row, with x- and y-coordinates) %%
%% Edges: list of neighbors (one vertex per row, at most six neighbors)    %%
%% Pairs: pairs of connected vertices (one edge per row)                   %%


dimns = dlmread(strcat(datapath,'.dimns'));
xmin = dimns(1,1);
xmax = dimns(1,2);
ymin = dimns(2,1);
ymax = dimns(2,2);
nDemes = xDemes*yDemes;
nPairs = nDemes*(nDemes-1)/2;
nEdges = (xDemes-1)*yDemes+(2*xDemes-1)*(yDemes-1);

%% A triangular grid extends half a triangle on the right %%
scalex = 1;
scaley = 1;
if ( (xmin<xmax) && (xDemes>1) ) 
  scalex = (xmax-xmin)/(xDemes-0.5);
end
if ( (ymin<ymax) && (yDemes>1) ) 
  scaley = (ymax-ymin)/(yDemes-1);
end

outer = [xmin,ymin;xmin,ymax;xmax,ymax;xmax,ymin;xmin,ymin];
habitat = struct('xmin',{xmin},'xmax',{xmax},...
                 'ymin',{ymin},'ymax',{ymax},...
		 'xDemes',{xDemes},...
		 'yDemes',{yDemes},...
		 'outer',{outer});
Demes = zeros(nDemes,2);
Edges = zeros(nDemes,6);
Pairs = zeros(nEdges*2,2);
e = 0;

for r = 1:yDemes
for c = 1:xDemes
  alpha = (r-1)*xDemes+c;
  Demes(alpha,1) = xmin+scalex*(c-1+0.5*mod(r-1,2));
  Demes(alpha,2) = ymin+scaley*(r-1);
  for pos = 1:6
    beta = fixed_neighbor(r,c,pos,xDemes,yDemes);
    Edges(alpha,pos) = beta;
    if (beta>0)
      e = e+1;
      Pairs(e,1) = alpha;
      Pairs(e,2) = beta;
    end
  end
end
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function beta = fixed_neighbor(r,c,pos,xDemes,yDemes)
%% In a triangular grid, a vertex has (at most) six
%% six neighbors. This function returns the indices
%% of the neighbors

% It is simpler to suppose indices start from 0
r = r-1;
c = c-1;
alpha = r*xDemes+c;
beta = -1;

if     ( (pos==1) && (mod(alpha  ,xDemes)>0) )
  beta = alpha-1;
elseif ( (pos==4) && (mod(alpha+1,xDemes)>0) )
  beta = alpha+1; 
elseif ( (pos==2) && (r<(yDemes-1) && mod(alpha  ,2*xDemes)>0) )
  beta = xDemes*(r+1)+c-mod(r+1,2);
elseif ( (pos==3) && (r<(yDemes-1) && mod(alpha+1,2*xDemes)>0) )
  beta = xDemes*(r+1)+c+1-mod(r+1,2);
elseif ( (pos==5) && (r>0 && mod(alpha+1,2*xDemes)>0) )
  beta = xDemes*(r-1)+c+1-mod(r+1,2);
elseif ( (pos==6) && (r>0 && mod(alpha  ,2*xDemes)>0) )
  beta = xDemes*(r-1)+c-  mod(r+1,2);
end

% In MATLAB indices start from 1
alpha = alpha+1;
beta = beta+1;
