

function [Graph,Data,Jindex] = suff_diploid_SNPs(datapath,Demes,Graph)
%% Read the sufficient statistic Diffs (the matrix of observed pairwise
%% genetic differences) and precompute some auxiliary matrices used to
%% compute the log likelihood
%% Diffs is precomputed but if Z = (Z_{ij}) is the allele count matrix
%% for i = 1..n individuals at j = 1..p SNPs, then
%%   the similarities are Sims = Z*Z'/p
%%   the dissimilarities are Diffs = ov*s0' + s0*ov' - 2*Sims
%%   where s0 = diag(Sims) are the self-similarities
%%   and ov is an n-vector of ones


Testing = 0;

dimns = dlmread(strcat(datapath,'.dimns'));
Coord = dlmread(strcat(datapath,'.coord'));
Diffs = dlmread(strcat(datapath,'.diffs'));
nIndiv = dimns(3,1);
nSites = dimns(3,2);

if nrow(Diffs)~=nIndiv || ncol(Diffs)~=nIndiv
  error('The dissimilarity matrix is not n-by-n.')
end
if nrow(Coord)~=nIndiv || ncol(Coord)~=2
  error('The coordinates matrix is not n-by-2.')
end
if ~issymetric(Diffs)
   error('The dissimilarity matrix is not symmetric.')
end
if ~isdistmat(Diffs)
  error('The dissimilarity matrix is not a valid distance matrix.')
end

[oDemes,Jinvpt,Jindex] = samples_to_demes(Coord,Demes);
n = length(Jinvpt);
o = length(oDemes);

Jpt = sparse(1:n,Jinvpt,1);
JtD = Jpt'*Diffs;
JtDJ = JtD*Jpt;
Sizes = Jpt'*ones(n,1);

%% A basis for contrasts %%
L = [-ones(n-1,1),eye(n-1)];
LDLt = - L*Diffs*L';
ldLLt = logdet(L*L');
ldLDLt = pseudologdet(LDLt);
ldDviQ = ldLLt - ldLDLt;

Data = struct('nIndiv',{nIndiv},...
              'nSites',{nSites},...
	      'Diffs',{Diffs},...
	      'Sizes',{Sizes},...
	      'microsat',{0},...
	      'diploid',{1},...
	      'JtDJ',{JtDJ},...
	      'ldLLt',{ldLLt},...
	      'ldDviQ',{ldDviQ},...
	      'Testing',{Testing});

if (Testing)
  Data.L = L;
  Data.oDemes = oDemes;
  Data.Jinvpt = Jinvpt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nDemes = nrow(Demes);
Juniq = oDemes;
Kuniq = 1:nDemes;
allobsrv = 1;
if length(Juniq)~=length(Kuniq)
  Kuniq(Juniq) = [];
  allobsrv = 0;
end
Graph.Juniq = Juniq;
Graph.Kuniq = Kuniq;
Graph.allobsrv = allobsrv;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
