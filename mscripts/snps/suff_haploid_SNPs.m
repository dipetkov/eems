

function [Mstruct,Sstruct,Jindex] = suff_haploid_SNPs(datapath,Demes,Mstruct)
%% Read the sufficient statistic Diffs (the matrix of observed pairwise %%
%% genetic differences) and precompute some auxiliary matrices used to  %%
%% compute the log likelihood                                           %%
%% Diffs is precomputed but if Z = (Z_{ij}) is the allele count matrix  %%
%% for i = 1..n individuals at j = 1..p SNPs, then                      %%
%%   the similarities are Sims = Z*Z'/p                                 %%
%%   the dissimilarities are Diffs = v*D' + D*v' - 2*Sims               %%
%%   where D = diag(Sims) and v is a vector of ones                     %%


%% ZZt is the _average_ similarity matrix
Coord = dlmread(strcat(datapath,'.coord'));
dimns = dlmread(strcat(datapath,'.dimns'));
Diffs = dlmread(strcat(datapath,'.diffs'));
nIndiv = dimns(3,1);
nSites = dimns(3,2);

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
winv = ones(o,1);

Jpt = sparse(1:n,Jinvpt,1);
JtO = Jpt'*ones(n);
JtD = Jpt'*Diffs;
JtOJ = JtO*Jpt;
JtDJ = JtD*Jpt;
Sizes = Jpt'*ones(n,1);
cwinvt = Sizes*winv';
JtDJvct = JtDJ*cwinvt' + cwinvt*JtDJ';

Bconst = winv'*JtDJ*winv;    %% winv'*J'*D*J*winv
oDinvoconst = - Sizes'*winv; %% -diag(J*J')'*winv
%% log(1*1') - logdet(J*J') + logdet(J*winv)
ldDinvconst = log(n) - sum(log(Sizes)) + sum(log(winv).*Sizes); 

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
                 'ldDinvconst',{ldDinvconst},...
                 'oDinvoconst',{oDinvoconst},...
                 'Bconst',{Bconst},...
		 'JtOJ',{JtOJ},...
		 'JtDJ',{JtDJ},...
		 'ldDviQ',{ldDviQ},...
		 'JtDJvct',{JtDJvct});
