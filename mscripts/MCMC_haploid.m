

function MCMC_haploid(sourcepath,datapath,mcmcpath,xPop,yPop,varargin)


addpath(strcat(sourcepath,'/mscripts/general'),...
        strcat(sourcepath,'/mscripts/project'),...
        strcat(sourcepath,'/mscripts/common'),...
        strcat(sourcepath,'/mscripts/mcmc'),...
        strcat(sourcepath,'/mscripts/snps'));

opt = default_params(sourcepath);

%%%%%%%% Process the input arguments %%%%%%%%
i = 0;
while (i<length(varargin))
  i = i + 1;
  if (strcmpi(varargin{i},'numMCMCIter'))
    opt.numMCMCIter = varargin{i+1};
  elseif (strcmpi(varargin{i},'numBurnIter'))
    opt.numBurnIter = varargin{i+1};
  elseif (strcmpi(varargin{i},'numThinIter'))
    opt.numThinIter = varargin{i+1};
  elseif (strcmpi(varargin{i},'mrateShape'))
    opt.mrateShape = varargin{i+1};
  elseif (strcmpi(varargin{i},'mrateScale'))
    opt.mrateScale = varargin{i+1};
  elseif (strcmpi(varargin{i},'s2locShape'))
    opt.s2locShape = varargin{i+1};
  elseif (strcmpi(varargin{i},'s2locScale'))
    opt.s2locScale = varargin{i+1};
  elseif (strcmpi(varargin{i},'negBiSize'))
    opt.negBiSize = varargin{i+1};
  elseif (strcmpi(varargin{i},'negBiProb'))
    opt.negBiProb = varargin{i+1};
  elseif (strcmpi(varargin{i},'mEffctProposalS2'))
    opt.mEffctProposalS2 = varargin{i+1};
  elseif (strcmpi(varargin{i},'mSeedsProposalS2'))
    opt.mSeedsProposalS2 = varargin{i+1};
  elseif (strcmpi(varargin{i},'mrateMuProposalS2'))
    opt.mrateMuProposalS2 = varargin{i+1};
  elseif (strcmpi(varargin{i},'dfProposalS2'))
    opt.dfProposalS2 = varargin{i+1};
  else
    error(['Invalid call to MCMC_microsats. The correct usage is:\n\t',...
	   'MCMC_haploid(sourcepath,datapath,mcmcpath,xPop,yPop,...)']);
  end
  i = i + 1;
end

[Demes,Edges,habitat,inPops,Mij] = make_triangular_grid(datapath,xPop,yPop);

%% Possibly remove vertex alpha from the graph by setting inPops(alpha) = 0 %%
%% The population graph must remain connected %%

[Mstruct,Demes,Edges] = make_irregular_grid(datapath,inPops,Mij,Demes,Edges);
[Mstruct,Sstruct,Jindex] = suff_haploid_SNPs(datapath,Demes,Mstruct);

gridsize = strcat(num2str(xPop),'x',num2str(yPop));
fprintf(2,'\nProcessing dataset %s\n',datapath);
fprintf(2,'MCMC output saved to %s\n',mcmcpath);
fprintf(2,'The triangular grid is %s\n',gridsize);
fprintf(2,'Input parameter values:\n');
dispstruct(2,opt);
fprintf(2,'\n\n');
  
[mVoronoi,params] = initial_values(Sstruct,Demes,habitat,opt);
mRates = realpow(10,params.mrateMu + mVoronoi.mEffcts);
mValues = average_rates(Mstruct,mRates,mVoronoi.mSeeds,mVoronoi.Demes);
kernel = resistance_kernel(Sstruct,Mstruct,mValues);

params = adjust_params(Sstruct,kernel,params);
initpi = wishart_prior(Sstruct,mVoronoi,params);
initll = wishart_ll(Sstruct,kernel,params);
[mcmc,schedule] = MCMC_initialize(mVoronoi,params,opt);

numSteps = mcmc.numSteps;
mcmcmhyper = zeros(numSteps,2);
mcmcthetas = zeros(numSteps,2);
mcmcpilogl = zeros(numSteps,2);
mcmcmtiles = zeros(numSteps,1);
mcmcmEffcts = cell(numSteps,1);
mcmcmRates = cell(numSteps,1);
mcmcxCoord = cell(numSteps,1);
mcmcyCoord = cell(numSteps,1);

fprintf(2,'Initial log prior:       %7.5f\n',initpi);
fprintf(2,'Initial log likelihood:  %7.5f\n',initll);
fprintf(2,'Initial log posterior:   %7.5f\n',initll+initpi);

pi0 = initpi;
ll0 = initll;

while ~mcmc.isfinished

  mcmc = MCMC_start_iteration(mcmc);

  while ~mcmc.iterDone

    %% Get a new configuration %%
    if (mcmc.currType==1)
      %% Propose to update the scale parameter sigma^22 and the degrees of freedom %%
      [proposal,pi1_pi0] = ...
          propose_thetas(Sstruct,kernel,params,schedule);
      if isinf(pi1_pi0) && (pi1_pi0<0)
        ll1 = 0;
      else
        ll1 = wishart_ll(Sstruct,kernel,proposal.params);
      end
    elseif (mcmc.currType==2)
      %% Propose to update the migration log rates %%
      [proposal,pi1_pi0] = ...
          propose_mEffects(kernel,params,mVoronoi,Sstruct,Mstruct,schedule);
      if isinf(pi1_pi0) && (pi1_pi0<0)
        ll1 = 0;
      else
        ll1 = wishart_ll(Sstruct,proposal.kernel,params);
      end
    elseif (mcmc.currType==3)
      %% Propose to update the mean migration log rate %%
      [proposal,pi1_pi0] = ...
          propose_mrateMu(kernel,params,mVoronoi,Sstruct,Mstruct);
      if isinf(pi1_pi0) && (pi1_pi0<0)
        ll1 = 0;
      else
        ll1 = wishart_ll(Sstruct,proposal.kernel,params);
      end
    elseif (mcmc.currType==4)
      %% Propose to move the Voronoi tiles around %%
      [proposal,pi1_pi0] = ...
          move_mVoronoi(kernel,params,mVoronoi,Sstruct,Mstruct,schedule);
      if isinf(pi1_pi0) && (pi1_pi0<0)
        ll1 = 0;
      else
        ll1 = wishart_ll(Sstruct,proposal.kernel,params);
      end
    elseif (mcmc.currType==5)
      %% Propose the birth or the death of a tile %%
      [proposal,pi1_pi0] = ...
          birthdeath_mVoronoi(kernel,params,mVoronoi,Sstruct,Mstruct);
      if isinf(pi1_pi0) && (pi1_pi0<0)
        ll1 = 0;
      else
        ll1 = wishart_ll(Sstruct,proposal.kernel,params);
      end
    end

    %% Metropolis-Hastings acceptance probability %%
    mcmc = MCMC_proceed(mcmc,proposal);
    u = rand(1);
    if isinf(pi1_pi0) && (pi1_pi0<0)
      %% Proposal is not accepted because the prior on it is 0 %%
      alpha = -Inf;
    else
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% Acceptance ratio (on the log scale)%%
      alpha = min(0,pi1_pi0+(ll1-ll0));     %%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (log(u)<alpha)
        [mVoronoi,kernel,params] = ...
            accept_move(mVoronoi,kernel,params,proposal);
        mcmc = MCMC_accept_proposal(mcmc); ll0 = ll1;
      end
    end

    if MCMC_exceeds_moves(mcmc,schedule)
      [mcmc,schedule] = index_schedule(mVoronoi,params,mcmc,schedule);
    end
    
    [mcmc,schedule] = MCMC_change_update(mcmc,schedule);

  end
  
  %% Update the hyperparameters (the cell effect variance tau_m^2)
  %% This changes the prior (but not the log likelihood)
  currmEffcts = mVoronoi.mEffcts;
  currmRates = realpow(10,params.mrateMu + currmEffcts);
  params = update_hyperp(Sstruct,mVoronoi,params,mcmc);
  params = update_dfsupp(Sstruct,mVoronoi,params,mcmc);
  pi0 = wishart_prior(Sstruct,mVoronoi,params);

  if sum(mcmc.currIter==mcmc.ii)
    currIndex = mcmc.kk(mcmc.currIter==mcmc.ii);
    mcmcmhyper(currIndex,:) = [params.mrateMu,params.mrateS2];
    mcmcthetas(currIndex,:) = [params.s2loc,params.df];
    mcmcpilogl(currIndex,:) = [pi0,ll0];
    mcmcmtiles(currIndex) = mVoronoi.mtiles;
    mcmcmEffcts{currIndex} = currmEffcts;
    mcmcmRates{currIndex} = currmRates;
    mcmcxCoord{currIndex} = mVoronoi.mSeeds(:,1);
    mcmcyCoord{currIndex} = mVoronoi.mSeeds(:,2);
  end

  mcmc = MCMC_end_iteration(mcmc,ll0);

  if mod(mcmc.currIter,100)==1
    fprintf(2,'    effective degrees of freedom = %.2f\n',params.df);
    fprintf(2,'        number of mVoronoi tiles = %d\n',mVoronoi.mtiles);
  end

end

save(strcat(mcmcpath,'.mat'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mcmc = MCMC_fdisp(2,mcmc);
fprintf(2,'Initial log prior:       %7.5f\n',initpi);
fprintf(2,'Initial log likelihood:  %7.5f\n',initll);
fprintf(2,'Final log prior:         %7.5f\n',pi0);
fprintf(2,'Final log likelihood:    %7.5f\n',ll0);

runid = fopen(strcat(mcmcpath,'.txt'),'w');
fprintf(runid,'Input parameter values:\n');
dispstruct(runid,opt);
fprintf(runid,'\nProcessing dataset %s\n',datapath);
fprintf(runid,'MCMC output saved to %s\n',mcmcpath);
fprintf(runid,'The triangular grid is %s\n',gridsize);
fprintf(runid,'Initial log prior:       %7.5f\n',initpi);
fprintf(runid,'Initial log likelihood:  %7.5f\n',initll);
fprintf(runid,'Final log prior:         %7.5f\n',pi0);
fprintf(runid,'Final log likelihood:    %7.5f\n',ll0);
mcmc = MCMC_fdisp(runid,mcmc);
fclose(runid);

dlmwrite(strcat(mcmcpath,'.mcmcmhyper'),mcmcmhyper,'delimiter',' ','precision',6);
dlmwrite(strcat(mcmcpath,'.mcmcthetas'),mcmcthetas,'delimiter',' ','precision',6);
dlmwrite(strcat(mcmcpath,'.mcmcpilogl'),mcmcpilogl,'delimiter',' ','precision',6);

dlmwrite(strcat(mcmcpath,'.mcmcmtiles'),mcmcmtiles,'delimiter',' ');

dlmcell(strcat(mcmcpath,'.mcmcxcoord'),mcmcxCoord,'delimiter',' ','precision',6);
dlmcell(strcat(mcmcpath,'.mcmcycoord'),mcmcyCoord,'delimiter',' ','precision',6);
dlmcell(strcat(mcmcpath,'.mcmcmrates'),mcmcmRates,'delimiter',' ','precision',6);

dlmwrite(strcat(mcmcpath,'.edges'),Edges,'delimiter',' ');
dlmwrite(strcat(mcmcpath,'.demes'),Demes,'delimiter',' ');
dlmwrite(strcat(mcmcpath,'.ipmap'),Jindex,'delimiter',' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fitted and observed distance matrices (for scatterplots)

dimns = dlmread(strcat(datapath,'.dimns'));
Coord = dlmread(strcat(datapath,'.coord'));
nIndiv = dimns(3,1);
nSites = dimns(3,2);
nDemes = nrow(Demes);

[oDemes,Jinvpt] = samples_to_demes(Coord,Demes);
Jpt = sparse(1:nIndiv,Jinvpt,1);
onev = ones(nIndiv,1);
Sizes = Jpt'*onev;
Counts = Sizes*Sizes';
Counts = Counts-diag(Sizes);

Dobs = Sstruct.Diffs;
Dhat = zeros(nIndiv);

nIters = length(mcmcmtiles);

for iter = 1:nIters

  mtiles = mcmcmtiles(iter);
  mRates = mcmcmRates{iter};
  xCoord = mcmcxCoord{iter};
  yCoord = mcmcyCoord{iter};
  mSeeds = [xCoord,yCoord];
  s2loc = mcmcthetas(iter,1);

  mValues = average_rates(Mstruct,mRates,mSeeds,Demes);
  qValues = ones(nDemes);

  G = sparse(Mstruct.Mi,Mstruct.Mj,mValues);
  R = resistance_distance(G);
  B = R(oDemes,oDemes)/4;
  w = qValues(oDemes);
  w = reshape(w,[],1);

  JBJt = Jpt*B*Jpt';
  Jw = Jpt*w;
  D = JBJt + onev*Jw'/2 + Jw*onev'/2 - diag(Jw);
  Dhat = Dhat + D*s2loc;

end

Dhat = Dhat/nIters;
JtDobsJ = Jpt'*Dobs*Jpt;
JtDhatJ = Jpt'*Dhat*Jpt;

%% Add 1 to the 0s to avoid dividing by 0
Counts = Counts + ~Counts;
JtDobsJ = JtDobsJ./Counts;
JtDhatJ = JtDhatJ./Counts;

dlmwrite(strcat(mcmcpath,'.rdistJtDobsJ'),JtDobsJ,'delimiter',' ','precision',6);
dlmwrite(strcat(mcmcpath,'.rdistJtDhatJ'),JtDhatJ,'delimiter',' ','precision',6);
