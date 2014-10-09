

function [Graph,Data,Jindex] = suff_microsats(datapath,Demes,Graph)
%% Read in the raw count data Z = (Z_{ik}) for i = 1..n individuals and %%
%% j = 1..p microsatellites and precompute some auxiliary matrices used %%
%% to compute the log likelihood                                        %%
%% There are two copies for each site, so Z has 2*p columns             %%
%% and missing alleles are coded as a negative number                   %%


Testing = 0;

%% Data is stored with one individual per line %%
Sites = dlmread(strcat(datapath,'.sites'));    %% alleles
Coord = dlmread(strcat(datapath,'.coord'));    %% coordinates
nSites = ncol(Sites)/2;

%% Pattern of missingness most likely varies across microsats
%% so save information for each locus in a cell array
%% Then will loop over the loci to compute the log likelihood (ll),
%% so wishart_ll takes longer to evaluate for microsats than for snps
Sizes = cell(nSites,1);
JtOJ = cell(nSites,1); 
JtDJ = cell(nSites,1);
JtDJvct = cell(nSites,1);
oDemes = cell(nSites,1);
nIndiv = zeros(nSites,1);
Bconst = zeros(nSites,1);
oDinvoconst = zeros(nSites,1);
ldDinvconst = zeros(nSites,1);

Z = cell(nSites,1);
O = cell(nSites,1);

for s = 1:nSites

  %% There are two alleles
  a1 = Sites(:,2*s-1);
  a2 = Sites(:,2*s);
  %% It shouldn't be the case that one allele
  %% is missing but the other is not
  if sum((a1<0)~=(a2<0))
    error('Error: one allele is missing but the other is not.')
  else
    obsrv = ~(a1<0);
  end

  n = sum(obsrv);
  a1 = a1(obsrv);
  a2 = a2(obsrv);
  z = (a1+a2)/2;
  nIndiv(s) = n;

  %% Similarity matrix (with rank 1)
  Sim = z*z';
  %% Dissimilarity matrix
  s0 = diag(Sim)*ones(1,n);
  Diff = s0 + s0' - 2*Sim;

  %% Demes: deme locations
  %% Samples: sampling locations
  Samples = Coord(obsrv,:);
  %% Assign each sample to the closest deme
  %% oDemes: indices of demes with at least one observation assigned
  %% oDemes(Jinvpt): indicates which deme each sample is assigned to
  [oDemes{s},Jinvpt] = samples_to_demes(Samples,Demes);

  n = length(Jinvpt);
  o = length(oDemes{s});

  %% Jpt is an indicator variable such that
  %% Jpt(i,a) = 1 if individual i comes from deme a and 0 otherwise
  Jpt = sparse(1:n,Jinvpt,1);
  JtO = Jpt'*ones(n);
  JtD = Jpt'*Diff;
  JtOJ{s} = JtO*Jpt;
  JtDJ{s} = JtD*Jpt;
  Sizes{s} = Jpt'*ones(n,1);
  
  Z{s} = z;
  O{s} = obsrv;
end

Data = struct('nIndiv',{nIndiv},...
              'nSites',{nSites},...
	      'oDemes',{oDemes},...
	      'Sizes',{Sizes},...
	      'microsat',{1},...
	      'JtOJ',{JtOJ},...
              'JtDJ',{JtDJ},...
	      'Testing',{Testing});

[oDemes,Jinvpt] = samples_to_demes(Coord,Demes);
Jindex = oDemes(Jinvpt);

if (Testing)
  Data.oDemes2 = oDemes;
  Data.Jinvpt2 = Jinvpt;
  Data.Z = Z;
  Data.O = O;
end
