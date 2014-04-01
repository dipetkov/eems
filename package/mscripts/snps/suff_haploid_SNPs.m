

function [Mstruct,Sstruct,Jindex] = suff_haploid_SNPs(datapath,Vcoord,Mstruct)
%% Read the sufficient statistic Diffs (the matrix of observed pairwise %%
%% genetic differences) and precompute some auxiliary matrices used to  %%
%% compute the log likelihood                                           %%
%% Diffs is precomputed but if Z = (Z_{ij}) is the allele count matrix  %%
%% for i = 1..n individuals at j = 1..p SNPs, then                      %%
%%   the similarities are Sims = Z*Z'/p                                 %%
%%   the dissimilarities are Diffs = v*D' + D*v' - 2*Sims               %%
%%   where D = diag(Sims) and v is a vector of ones                     %%


%% Important:
%% ZZt is the _average_ similarity matrix
Jcoord = dlmread(strcat(datapath,'.coord'));
dimns = dlmread(strcat(datapath,'.dimns'));
Diffs = dlmread(strcat(datapath,'.diffs'));
nIndiv = dimns(3,1);
nSites = dimns(3,2);

if rank(Diffs)<nIndiv
  error('The dissimilarity matrix is rank-deficient.')
end

[tempi,Jindex] = pdist2(Vcoord,Jcoord,'euclidean');
[Juniq,tempi,Jinvpt] = unique(Jindex);
[dC,u] = hist(Jinvpt,unique(Jinvpt));
nDemes = Mstruct.nDemes;
oDemes = length(Juniq);

Kuniq = 1:nDemes;
allobsrv = 1;
if (oDemes~=nDemes)
  Kuniq(Juniq) = [];
  allobsrv = 0;
end

Mstruct.Juniq = Juniq;
Mstruct.Kuniq = Kuniq;
Mstruct.allobsrv = allobsrv;

n = nIndiv;
o = oDemes;
oDinvoconst = - n;
ldDinvconst = + log(n);
Bconst = sum(sum(Diffs));

Onev = ones(n);
Jpt = sparse(1:n,Jinvpt,1);
JtOJ = Jpt'*Onev*Jpt;
JtDJ = Jpt'*Diffs*Jpt;
JtDOJ = Jpt'*(Diffs*Onev+Onev*Diffs)*Jpt;
Counts = sparse(1:o,1:o,dC);
Cinv = sparse(1:o,1:o,1./dC);

DviQ = projection_DinvQ(Diffs/2);
ldDviQ = pseudo_logdet(-DviQ);

Sstruct = struct('nIndiv',{nIndiv},...
                 'nSites',{nSites},...
                 'oDemes',{oDemes},...
                 'Jinvpt',{Jinvpt},...
                 'Counts',{Counts},...
                 'ldDinvconst',{ldDinvconst},...
                 'oDinvoconst',{oDinvoconst},...
                 'Bconst',{Bconst},...
		 'microsat',{0},...
                 'diploid',{0},...
		 'Cconst',{1},...
		 'Cinv',{Cinv},...
		 'JtOJ',{JtOJ},...
		 'JtDJ',{JtDJ},...
		 'JtDOJ',{JtDOJ},...
		 'ldDviQ',{ldDviQ});
