

function MCMC_haploid(sourcepath,datapath,mcmcpath,xPop,yPop,varargin)


addpath(strcat(sourcepath,'/mscripts/general'),...
        strcat(sourcepath,'/mscripts/project'),...
        strcat(sourcepath,'/mscripts/common'),...
        strcat(sourcepath,'/mscripts/mcmc'),...
        strcat(sourcepath,'/mscripts/snps'),...
        option='-end');

opt = default_snps_params(sourcepath);

%%%%%%%% Process the input arguments %%%%%%%%
i = 0;
while (i<length(varargin))
  i++;
  if (strcmpi(varargin{i},"numIter"))
    opt.numIter = varargin{++i};
  elseif (strcmpi(varargin{i},"numBurn"))
    opt.numBurn = varargin{++i};
  elseif (strcmpi(varargin{i},"numThin"))
    opt.numThin = varargin{++i};
  elseif (strcmpi(varargin{i},"effctS2"))
    opt.effctS2 = varargin{++i};
  elseif (strcmpi(varargin{i},"coordS2"))
    opt.coordS2 = varargin{++i};
  elseif (strcmpi(varargin{i},"log10S2"))
    opt.log10S2 = varargin{++i};
  elseif (strcmpi(varargin{i},"ratesC"))
    opt.ratesC = varargin{++i};
  elseif (strcmpi(varargin{i},"s2locC"))
    opt.s2locC = varargin{++i};
  elseif (strcmpi(varargin{i},"negBiR"))
    opt.negBiR = varargin{++i};
  elseif (strcmpi(varargin{i},"negBiP"))
    opt.negBiP = varargin{++i};
  elseif (strcmpi(varargin{i},"dfS2"))
    opt.dfS2 = varargin{++i};
  else
    print_usage( );
  end
end

[Vcoord,Vedges,Mij] = make_triangular_grid(datapath,xPop,yPop);
[habitat,inPops] = read_rectangular_habitat(datapath,Vcoord,xPop,yPop);

%% Possibly remove vertex alpha from the graph by setting inPops(alpha) = 0 %%
%% The population graph must remain connected %%

[Mstruct,Vcoord,Vedges] = make_irregular_grid(datapath,inPops,Mij,Vcoord,Vedges);
[Mstruct,Sstruct,Jindex] = suff_haploid_SNPs(datapath,Vcoord,Mstruct);

gridsize = strcat(num2str(xPop),'x',num2str(yPop));
fprintf(2,'\nProcessing dataset %s\n',datapath);
fprintf(2,'MCMC output saved to %s\n',mcmcpath);
fprintf(2,'The triangular grid is %s\n',gridsize);
fprintf(2,'Input parameter values:\n');
fdisp(2,opt);
fprintf(2,'\n\n');
  
[Voronoi,params] = assign_effects(Sstruct,Vcoord,habitat,opt);
mRates = realpow(10,params.rateMu + Voronoi.Effcts);
mValues = average_rates(Mstruct,mRates,Voronoi.Scoord,Voronoi.Vcoord);
kernel = resistance_kernel(Sstruct,Mstruct,mValues);

params = adjust_params(Sstruct,kernel,params);
initpi = wishart_prior(Sstruct,Voronoi,params);
[initll,inittr] = wishart_ll(Sstruct,kernel,params);
[mcmc,schedule] = MCMC_initialize(Voronoi,params,opt);
params.trDinvQxD = inittr;

numSteps = mcmc.numSteps;
mcmchyperP = zeros(numSteps,2);
mcmcthetas = zeros(numSteps,2);
mcmcpilogl = zeros(numSteps,2);
mcmcNtiles = zeros(numSteps,1);
mcmcEffcts = cell(numSteps,1);
mcmcmRates = cell(numSteps,1);
mcmcXcoord = cell(numSteps,1);
mcmcYcoord = cell(numSteps,1);

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
        [ll1,tr1] = wishart_ll(Sstruct,kernel,proposal.params);
      end
    elseif (mcmc.currType==2)
      %% Propose to update the migration log rates %%
      [proposal,pi1_pi0] = ...
          propose_effects(kernel,params,Voronoi,Sstruct,Mstruct,schedule);
      if isinf(pi1_pi0) && (pi1_pi0<0)
        ll1 = 0;
      else
        [ll1,tr1] = wishart_ll(Sstruct,proposal,params);
      end
    elseif (mcmc.currType==3)
      %% Propose to update the mean migration log rate %%
      [proposal,pi1_pi0] = ...
          propose_rateMu(kernel,params,Voronoi,Sstruct,Mstruct);
      if isinf(pi1_pi0) && (pi1_pi0<0)
        ll1 = 0;
      else
        [ll1,tr1] = wishart_ll(Sstruct,proposal,params);
      end
    elseif (mcmc.currType==4)
      %% Propose to move the Voronoi tiles around %%
      [proposal,pi1_pi0] = ...
          move_nuclei(kernel,params,Voronoi,Sstruct,Mstruct,schedule);
      if isinf(pi1_pi0) && (pi1_pi0<0)
        ll1 = 0;
      else
        [ll1,tr1] = wishart_ll(Sstruct,proposal,params);
      end
    elseif (mcmc.currType==5)
      %% Propose the birth or the death of a tile %%
      [proposal,pi1_pi0] = ...
          propose_birth_death(kernel,params,Voronoi,Sstruct,Mstruct);
      if isinf(pi1_pi0) && (pi1_pi0<0)
        ll1 = 0;
      else
        [ll1,tr1] = wishart_ll(Sstruct,proposal,params);
      end
    end

    %% Metropolis-Hastings acceptance probability %%
    mcmc = MCMC_proceed(mcmc,proposal);
    u = unifrnd(0,1);
    if isinf(pi1_pi0) && (pi1_pi0<0)
      %% Proposal is not accepted because the prior on it is 0 %%
      alpha = -Inf;
    else
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% Acceptance ratio (on the log scale)%%
      alpha = min(0,pi1_pi0+(ll1-ll0));     %%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (log(u)<alpha)
        [Voronoi,kernel,params] = ...
            accept_move(Voronoi,kernel,params,proposal);
        ll0 = ll1; params.trDinvQxD = tr1;
        mcmc = MCMC_accept_proposal(mcmc);
      end
    end

    if MCMC_exceeds_moves(mcmc,schedule)
      [mcmc,schedule] = index_schedule(Voronoi,params,mcmc,schedule);
    end
    
    [mcmc,schedule] = MCMC_change_update(mcmc,schedule);

  end
  
  %% Update the hyperparameters (the cell effect variance tau_m^2)
  %% This changes the prior (but not the log likelihood)
  currEffcts = Voronoi.Effcts;
  currmRates = realpow(10,params.rateMu + currEffcts);
  params = update_hyperp(Sstruct,Voronoi,params,mcmc);
  params = update_dfsupp(Sstruct,Voronoi,params,mcmc);
  pi0 = wishart_prior(Sstruct,Voronoi,params);

  if sum(mcmc.currIter==mcmc.ii)
    currIndex = mcmc.kk(mcmc.currIter==mcmc.ii);
    mcmchyperP(currIndex,:) = [params.rateMu,params.rateS2];
    mcmcthetas(currIndex,:) = [params.s2loc,params.df];
    mcmcpilogl(currIndex,:) = [pi0,ll0];
    mcmcNtiles(currIndex) = Voronoi.ntiles;
    mcmcEffcts{currIndex} = currEffcts;
    mcmcmRates{currIndex} = currmRates;
    mcmcXcoord{currIndex} = Voronoi.Scoord(:,1);
    mcmcYcoord{currIndex} = Voronoi.Scoord(:,2);
  end

  mcmc = MCMC_end_iteration(mcmc,ll0);

  if mod(mcmc.currIter,100)==1
    fprintf(2,'    effective degrees of freedom = %.2f\n',params.df);
    fprintf(2,'         number of Voronoi tiles = %d\n',Voronoi.ntiles);
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
fdisp(runid,opt);
fprintf(runid,'\nProcessing dataset %s\n',datapath);
fprintf(runid,'MCMC output saved to %s\n',mcmcpath);
fprintf(runid,'The triangular grid is %s\n',gridsize);
fprintf(runid,'Initial log prior:       %7.5f\n',initpi);
fprintf(runid,'Initial log likelihood:  %7.5f\n',initll);
fprintf(runid,'Final log prior:         %7.5f\n',pi0);
fprintf(runid,'Final log likelihood:    %7.5f\n',ll0);
mcmc = MCMC_fdisp(runid,mcmc);
fclose(runid);

dlmwrite(strcat(mcmcpath,'.mcmchyperp'),mcmchyperP,'delimiter',' ');
dlmwrite(strcat(mcmcpath,'.mcmcthetas'),mcmcthetas,'delimiter',' ');
dlmwrite(strcat(mcmcpath,'.mcmcntiles'),mcmcNtiles,'delimiter',' ');
dlmwrite(strcat(mcmcpath,'.mcmcpilogl'),mcmcpilogl,'delimiter',' ');

dlmcell(strcat(mcmcpath,'.mcmcxcoord'),mcmcXcoord);
dlmcell(strcat(mcmcpath,'.mcmcycoord'),mcmcYcoord);
dlmcell(strcat(mcmcpath,'.mcmceffcts'),mcmcEffcts);
dlmcell(strcat(mcmcpath,'.mcmcmrates'),mcmcmRates);

dlmwrite(strcat(mcmcpath,'.edges'),Vedges,'delimiter',' ');
dlmwrite(strcat(mcmcpath,'.coord'),Vcoord,'delimiter',' ');
dlmwrite(strcat(mcmcpath,'.ipmap'),Jindex,'delimiter',' ');
