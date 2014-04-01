

function [Mstruct,Sstruct,Jindex] = suff_microsats(datapath,Vcoord,Mstruct)
%% Read in the raw count data Z = (Z_{ik}) for i = 1..n individuals and %%
%% j = 1..p microsatellites and precompute some auxiliary matrices used %%
%% to compute the log likelihood                                        %%
%% There are two copies for each site, so Z has 2*p columns             %%
%% and missing alleles are coded as a negative number                   %%


%% Important:                                  %%
%% Data is stored with one individual per line %%
Zall = dlmread(strcat(datapath,'.sites'));
Jall = dlmread(strcat(datapath,'.coord'));

nSites = ncol(Zall)/2;
nIndiv = nrow(Zall);
Zstar = NaN(nIndiv,nSites);

Cinv = cell(nSites,1);
Counts = cell(nSites,1);
Bconst = zeros(nSites,1);
oDinvoconst = zeros(nSites,1);
ldDinvconst = zeros(nSites,1);
JtOJ = cell(nSites,1);
JtDJ = cell(nSites,1);
JtDOJ = cell(nSites,1);
nIndiv = zeros(nSites,1);
oDemes = zeros(nSites,1);
nDemes = Mstruct.nDemes;

Mstruct.Juniq = cell(nSites,1);
Mstruct.Kuniq = cell(nSites,1);
Mstruct.allobsrv = zeros(nSites,1);

for s = 1:nSites

  %% There are two haplotypes                 %%
  h1 = Zall(:,2*s-1);
  h2 = Zall(:,2*s);
  %% It shouldn't be the case that one allele %%
  %% is missing but the other is not          %%
  if sum((h1<0)~=(h2<0))
    error('Error: one allele is missing but the other is not.')
  end

  obsrv = (h1>0);
  n = sum(obsrv);
  h1 = h1(obsrv);
  h2 = h2(obsrv);

  nIndiv(s) = n;
  Jcoord = Jall(obsrv,:);
  [tempi,Jindex] = pdist2(Vcoord,Jcoord,'euclidean');
  [Juniq,tempi,Jinvpt] = unique(Jindex);
  [dC,u] = hist(Jinvpt,unique(Jinvpt));
  oDemes(s) = length(Juniq);

  Kuniq = 1:nDemes;
  allobsrv = 1;
  if (oDemes(s)~=nDemes)
    Kuniq(Juniq) = [];
    allobsrv = 0;
  end

  Mstruct.Juniq{s} = Juniq;
  Mstruct.Kuniq{s} = Kuniq;
  Mstruct.allobsrv(s) = allobsrv;

  zstar = (h1+h2)/2;
  Sims = zstar*zstar';
  %% Dissimilarity matrix %%
  dS = repmat(diag(Sims),1,n);
  Diffs = dS + dS' - 2*Sims;

  oDinvoconst(s) = - n;
  ldDinvconst(s) = + log(n);
  Bconst(s) = sum(sum(Diffs));

  o = oDemes(s);
  Onev = ones(n);
  Jpt = sparse(1:n,Jinvpt,1);
  JtOJ{s} = Jpt'*Onev*Jpt;
  JtDJ{s} = Jpt'*Diffs*Jpt;
  JtDOJ{s} = Jpt'*(Diffs*Onev+Onev*Diffs)*Jpt;
  Counts{s} = sparse(1:o,1:o,dC);
  Cinv{s} = sparse(1:o,1:o,1./dC);

  Zstar(obsrv,s) = zstar;
end

[temp,Jindex] = pdist2(Vcoord,Jall,'euclidean');

Sstruct = struct('nIndiv',{nIndiv},...
                 'nSites',{nSites},...
                 'oDemes',{oDemes},...
                 'diploid',{0},...
		 'microsat',{1},...
                 'Zstar',{Zstar},...
                 'ldDinvconst',{ldDinvconst},...
                 'oDinvoconst',{oDinvoconst},...
                 'Jinvpt',{Jinvpt},...
                 'Counts',{Counts},...
                 'Bconst',{Bconst},...
                 'Cinv',{Cinv},...
		 'JtOJ',{JtOJ},...
                 'JtDJ',{JtDJ},...
                 'JtDOJ',{JtDOJ});

