

function [Mstruct,Sstruct,Jindex] = suff_haploid_SNPs(datapath,Demes,Mstruct)
%% Read the sufficient statistic Diffs (the matrix of observed pairwise %%
%% genetic differences) and precompute some auxiliary matrices used to  %%
%% compute the log likelihood                                           %%
%% Diffs is precomputed but if Z = (Z_{ij}) is the allele count matrix  %%
%% for i = 1..n individuals at j = 1..p SNPs, then                      %%
%%   the similarities are Sims = Z*Z'/p                                 %%
%%   the dissimilarities are Diffs = ov*s0' + s0*ov' - 2*Sims           %%
%%   where s0 = diag(Sims) are the self-similarities                    %%
%%   and ov is an n-vector of ones                                      %%


Coord = dlmread(strcat(datapath,'.coord'));
dimns = dlmread(strcat(datapath,'.dimns'));
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
if rank(Diffs)<nIndiv
  error('The dissimilarity matrix is rank-deficient.')
end

[oDemes,Jinvpt] = samples_to_demes(Coord,Demes);
Jindex = oDemes(Jinvpt);

n = length(Jinvpt);
o = length(oDemes);

Jpt = sparse(1:n,Jinvpt,1);
JtO = Jpt'*ones(n);
JtD = Jpt'*Diffs;
JtOJ = JtO*Jpt;
JtDJ = JtD*Jpt;
Sizes = Jpt'*ones(n,1);

%% A basis for contrasts %%
L = [-ones(n-1,1),eye(n-1)];
ldDviQ = logdet(L*L') - logdet(-L*Diffs*L');

Sstruct = struct('nIndiv',{nIndiv},...
                 'nSites',{nSites},...
		 'oDemes',{oDemes},...
		 'Diffs',{Diffs},...
		 'Sizes',{Sizes},...
		 'microsat',{0},...
		 'diploid',{0},...
		 'JtOJ',{JtOJ},...
		 'JtDJ',{JtDJ},...
		 'ldDviQ',{ldDviQ});
