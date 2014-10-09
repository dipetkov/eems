

function MCMC_microsat(sourcepath,datapath,mcmcpath,xDemes,yDemes,varargin)


addpath(strcat(sourcepath,'/mscripts/general'),...
        strcat(sourcepath,'/mscripts/project'),...
        strcat(sourcepath,'/mscripts/common'),...
        strcat(sourcepath,'/mscripts/mcmc'),...
        strcat(sourcepath,'/mscripts/sats'));

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
  elseif (strcmpi(varargin{i},'qEffctProposalS2'))
    opt.qEffctProposalS2 = varargin{i+1};
  elseif (strcmpi(varargin{i},'mSeedsProposalS2'))
    opt.mSeedsProposalS2 = varargin{i+1};
  elseif (strcmpi(varargin{i},'qSeedsProposalS2'))
    opt.qSeedsProposalS2 = varargin{i+1};
  elseif (strcmpi(varargin{i},'mrateMuProposalS2'))
    opt.mrateMuProposalS2 = varargin{i+1};
  else
    error(['Invalid call to MCMC_microsat. The correct usage is:\n\t',...
	   'MCMC_microsat(sourcepath,datapath,mcmcpath,xDemes,yDemes,...)']);
  end
  i = i + 1;
end

[Demes,Edges,Pairs,habitat,inDemes] = make_triangular_grid(datapath,xDemes,yDemes);

%% Possibly remove vertex alpha from the graph by setting inDemes(alpha) = 0
%% The graph must remain connected

[Graph,Demes,Edges] = update_irregular_grid(datapath,inDemes,Demes,Edges,Pairs);
[Graph,Data,Jindex] = suff_microsats(datapath,Demes,Graph);

gridsize = strcat(num2str(xDemes),'x',num2str(yDemes));
fprintf(2,'\nProcessing dataset %s\n',datapath);
fprintf(2,'MCMC output saved to %s\n',mcmcpath);
fprintf(2,'The triangular grid is %s\n',gridsize);
fprintf(2,'Input parameter values:\n');
dispstruct(2,opt);
fprintf(2,'\n\n');
  
[qVoronoi,mVoronoi,params] = initial_values(Data,Demes,habitat,opt);
mRates = realpow(10,mVoronoi.mEffcts + params.mrateMu);
qRates = realpow(10,qVoronoi.qEffcts + params.qrateMu);
[qValues,mValues] = ...
  average_rates(Graph,qRates,mRates,qVoronoi.qSeeds,mVoronoi.mSeeds,Demes);
kernel = resistance_kernel(Data,Graph,mValues,qValues);

params = adjust_params(Data,kernel,params);
initpi = wishart_prior(Data,qVoronoi,mVoronoi,params);
initll = wishart_ll(Data,kernel,params);
[mcmc,schedule] = MCMC_initialize(qVoronoi,mVoronoi,params,opt);

numSteps = mcmc.numSteps;
mcmcmhyper = zeros(numSteps,2);
mcmcqhyper = zeros(numSteps,2);
mcmcthetas = zeros(numSteps,2);
mcmcpilogl = zeros(numSteps,2);
mcmcmtiles = zeros(numSteps,1);
mcmcqtiles = zeros(numSteps,1);
mcmcmEffcts = cell(numSteps,1);
mcmcqEffcts = cell(numSteps,1);
mcmcmRates = cell(numSteps,1);
mcmcqRates = cell(numSteps,1);
mcmcxCoord = cell(numSteps,1);
mcmcyCoord = cell(numSteps,1);
mcmcwCoord = cell(numSteps,1);
mcmczCoord = cell(numSteps,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mcmcvS2loc = zeros(numSteps,Data.nSites);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
      %% Propose to update the scale parameters sigma_s^2 %%
      %% There is one parameter for each microsatellite   %%
      [proposal,pi1_pi0] = ...
        propose_thetas(Data,kernel,params,schedule);
      if isinf(pi1_pi0) && (pi1_pi0<0)
        ll1 = 0;
      else
        ll1 = wishart_ll(Data,kernel,proposal.params);
      end
    elseif (mcmc.currType==2)
      %% Propose to update the coalescence log rates %%
      [proposal,pi1_pi0] = ...
        propose_qEffects(kernel,params,qVoronoi,mVoronoi,Data,Graph,schedule);
      if isinf(pi1_pi0) && (pi1_pi0<0)
        ll1 = 0;
      else
        ll1 = wishart_ll(Data,proposal.kernel,params);
      end
    elseif (mcmc.currType==3)
      %% Propose to update the mean coalescence log rate %%
      [proposal,pi1_pi0] = ...
        propose_qrateMu(kernel,params,qVoronoi,mVoronoi,Data,Graph);
      if isinf(pi1_pi0) && (pi1_pi0<0)
        ll1 = 0;
      else
        ll1 = wishart_ll(Data,proposal.kernel,params);
      end
    elseif (mcmc.currType==4)
      %% Propose to move the Voronoi tiles around %%
      [proposal,pi1_pi0] = ...
        move_qVoronoi(kernel,params,qVoronoi,mVoronoi,Data,Graph,schedule);
      if isinf(pi1_pi0) && (pi1_pi0<0)
        ll1 = 0;
      else
        ll1 = wishart_ll(Data,proposal.kernel,params);
      end
    elseif (mcmc.currType==5)
      %% Propose the birth or the death of a tile %%
      [proposal,pi1_pi0] = ...
        birthdeath_qVoronoi(kernel,params,qVoronoi,mVoronoi,Data,Graph);
      if isinf(pi1_pi0) && (pi1_pi0<0)
        ll1 = 0;
      else
        ll1 = wishart_ll(Data,proposal.kernel,params);
      end
    elseif (mcmc.currType==6)
      %% Propose to update the migration log rates %%
      [proposal,pi1_pi0] = ...
        propose_mEffects(kernel,params,qVoronoi,mVoronoi,Data,Graph,schedule);
      if isinf(pi1_pi0) && (pi1_pi0<0)
        ll1 = 0;
      else
        ll1 = wishart_ll(Data,proposal.kernel,params);
      end
    elseif (mcmc.currType==7)
      %% Propose to update the mean migration log rate %%
      [proposal,pi1_pi0] = ...
        propose_mrateMu(kernel,params,qVoronoi,mVoronoi,Data,Graph);
      if isinf(pi1_pi0) && (pi1_pi0<0)
        ll1 = 0;
      else
        ll1 = wishart_ll(Data,proposal.kernel,params);
      end
    elseif (mcmc.currType==8)
      %% Propose to move the Voronoi tiles around %%
      [proposal,pi1_pi0] = ...
        move_mVoronoi(kernel,params,qVoronoi,mVoronoi,Data,Graph,schedule);
      if isinf(pi1_pi0) && (pi1_pi0<0)
        ll1 = 0;
      else
        ll1 = wishart_ll(Data,proposal.kernel,params);
      end
    elseif (mcmc.currType==9)
      %% Propose the birth or the death of a tile %%
      [proposal,pi1_pi0] = ...
        birthdeath_mVoronoi(kernel,params,qVoronoi,mVoronoi,Data,Graph);
      if isinf(pi1_pi0) && (pi1_pi0<0)
        ll1 = 0;
      else
        ll1 = wishart_ll(Data,proposal.kernel,params);
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
        [qVoronoi,mVoronoi,kernel,params] = ...
            accept_move(qVoronoi,mVoronoi,kernel,params,proposal);
	mcmc = MCMC_accept_proposal(mcmc); ll0 = ll1;
      end
    end

    if MCMC_exceeds_moves(mcmc,schedule)
      [mcmc,schedule] = index_schedule(qVoronoi,mVoronoi,params,mcmc,schedule);
    end
    
    [mcmc,schedule] = MCMC_change_update(mcmc,schedule);

  end
  
  %% Update the hyperparameters (the effect variance)          %%
  %% This also changes the prior (but not the likelihood term) %%
  mRates = realpow(10,mVoronoi.mEffcts + params.mrateMu);
  qRates = realpow(10,qVoronoi.qEffcts + params.qrateMu);	
  params = update_hyperp(Data,qVoronoi,mVoronoi,params,mcmc);
  params = update_dfsupp(Data,params,mcmc);
  pi0 = wishart_prior(Data,qVoronoi,mVoronoi,params);

  if sum(mcmc.currIter==mcmc.ii)
    currIndex = mcmc.kk(mcmc.currIter==mcmc.ii);
    mcmcmhyper(currIndex,:) = [params.mrateMu,params.mrateS2];
    mcmcqhyper(currIndex,:) = [params.qrateMu,params.qrateS2];
    mcmcthetas(currIndex,:) = [-9,params.df];
    mcmcpilogl(currIndex,:) = [pi0,ll0];
    mcmcvS2loc(currIndex,:) = params.s2loc';
    mcmcmtiles(currIndex) = mVoronoi.mtiles;
    mcmcqtiles(currIndex) = qVoronoi.qtiles;
    mcmcmEffcts{currIndex} = mVoronoi.mEffcts;
    mcmcqEffcts{currIndex} = qVoronoi.qEffcts;
    mcmcmRates{currIndex} = mRates;
    mcmcqRates{currIndex} = qRates;
    mcmcxCoord{currIndex} = mVoronoi.mSeeds(:,1);
    mcmcyCoord{currIndex} = mVoronoi.mSeeds(:,2);
    mcmcwCoord{currIndex} = qVoronoi.qSeeds(:,1);
    mcmczCoord{currIndex} = qVoronoi.qSeeds(:,2);
  end

  mcmc = MCMC_end_iteration(mcmc,ll0);

  if mod(mcmc.currIter,100)==1
    fprintf(2,'    effective degrees of freedom = %.2f\n',params.df);
    fprintf(2,'        number of qVoronoi tiles = %d\n',qVoronoi.qtiles);
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

dlmwrite(strcat(mcmcpath,'.mcmcmhyper'),mcmcmhyper,'delimiter',' ','precision','%.6f');
dlmwrite(strcat(mcmcpath,'.mcmcqhyper'),mcmcqhyper,'delimiter',' ','precision','%.6f');
dlmwrite(strcat(mcmcpath,'.mcmcthetas'),mcmcthetas,'delimiter',' ','precision','%.6f');
dlmwrite(strcat(mcmcpath,'.mcmcpilogl'),mcmcpilogl,'delimiter',' ','precision','%.6f');
dlmwrite(strcat(mcmcpath,'.mcmcvS2loc'),mcmcvS2loc,'delimiter',' ','precision','%.6f');

dlmwrite(strcat(mcmcpath,'.mcmcmtiles'),mcmcmtiles,'delimiter',' '); %% integers
dlmwrite(strcat(mcmcpath,'.mcmcqtiles'),mcmcqtiles,'delimiter',' '); %% integers

dlmcell(strcat(mcmcpath,'.mcmcxcoord'),mcmcxCoord,'delimiter',' ','precision','%.6f');
dlmcell(strcat(mcmcpath,'.mcmcycoord'),mcmcyCoord,'delimiter',' ','precision','%.6f');
dlmcell(strcat(mcmcpath,'.mcmcwcoord'),mcmcwCoord,'delimiter',' ','precision','%.6f');
dlmcell(strcat(mcmcpath,'.mcmczcoord'),mcmczCoord,'delimiter',' ','precision','%.6f');
dlmcell(strcat(mcmcpath,'.mcmcmrates'),mcmcmRates,'delimiter',' ','precision','%.6f');
dlmcell(strcat(mcmcpath,'.mcmcqrates'),mcmcqRates,'delimiter',' ','precision','%.6f');

dlmwrite(strcat(mcmcpath,'.demes'),Demes,'delimiter',' ','precision','%.6f');
dlmwrite(strcat(mcmcpath,'.edges'),Edges,'delimiter',' '); %% integers
dlmwrite(strcat(mcmcpath,'.ipmap'),Jindex,'delimiter',' '); %% integers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fitted and observed distance matrices (for scatterplots)

Sites = dlmread(strcat(datapath,'.sites'));
Coord = dlmread(strcat(datapath,'.coord'));
nSites = ncol(Sites)/2;
nIndiv = nrow(Sites);
nDemes = nrow(Demes);

[oDemes,Jinvpt] = samples_to_demes(Coord,Demes);
Jpt = sparse(1:nIndiv,Jinvpt,1);
onev = ones(nIndiv,1);
Sizes = Jpt'*onev;
Counts = Sizes*Sizes';
Counts = Counts-diag(Sizes);

Dobs = zeros(nIndiv);
Dhat = zeros(nIndiv);
Pairs = zeros(nIndiv);

for s = 1:nSites
    
  a1 = Sites(:,2*s-1);
  a2 = Sites(:,2*s);
  obsrv = ~(a1<0)+0;
  a1(~obsrv) = NaN;
  a2(~obsrv) = NaN;

  z = (a1+a2)/2;
  z = z-nanmean(z);
  z = z/nanstd(z);
  z(isnan(z)) = 0;
  %% Similarity matrix (with rank 1)
  Sim = z*z';
  %% Dissimilarity matrix
  s0 = repmat(diag(Sim),1,nIndiv);
  Dobs = Dobs + s0 + s0' - 2*Sim;  
  Pairs = Pairs + obsrv*obsrv';

end

nIters = length(mcmcmtiles);

for iter = 1:nIters

  mRates = mcmcmRates{iter};
  qRates = mcmcqRates{iter};
  xCoord = mcmcxCoord{iter};
  yCoord = mcmcyCoord{iter};
  wCoord = mcmcwCoord{iter};
  zCoord = mcmczCoord{iter};
  mSeeds = [xCoord,yCoord];
  qSeeds = [wCoord,zCoord];

  [qValues,mValues] = average_rates(Graph,qRates,mRates,qSeeds,mSeeds,Demes);

  G = sparse(Graph.Vi,Graph.Vj,mValues);
  R = resistance_distance(G);
  B = R(oDemes,oDemes)/4; %% between-deme distances
  w = qValues(oDemes);    %% within-deme distances
  w = reshape(w,[],1);

  %% Construct the entire n-ny-n fitted distance matrix
  JBJt = Jpt*B*Jpt';
  Jw = Jpt*w;
  D = JBJt + onev*Jw'/2 + Jw*onev'/2 - diag(Jw);
  Dhat = Dhat + D;
  
end

Dobs = Dobs./Pairs;
Dhat = Dhat/nIters;
JtDobsJ = Jpt'*Dobs*Jpt;
JtDhatJ = Jpt'*Dhat*Jpt;

%% Add 1 to the 0s to avoid dividing by 0
Counts = Counts + ~Counts;
JtDobsJ = JtDobsJ./Counts;
JtDhatJ = JtDhatJ./Counts;

%dlmwrite(strcat(mcmcpath,'.rdistDobs'),Dobs,'delimiter',' ','precision','%.6f');
%dlmwrite(strcat(mcmcpath,'.rdistDhat'),Dhat,'delimiter',' ','precision','%.6f');
dlmwrite(strcat(mcmcpath,'.rdistSizes'),Sizes,'delimiter',' ');
dlmwrite(strcat(mcmcpath,'.rdistoDemes'),oDemes','delimiter',' ');
dlmwrite(strcat(mcmcpath,'.rdistJtDobsJ'),JtDobsJ,'delimiter',' ','precision','%.6f');
dlmwrite(strcat(mcmcpath,'.rdistJtDhatJ'),JtDhatJ,'delimiter',' ','precision','%.6f');
