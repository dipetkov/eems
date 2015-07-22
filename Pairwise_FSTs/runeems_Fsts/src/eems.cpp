
#include "eems.hpp"

EEMS::EEMS(const Params &params) {
  this->params = params;
  draw.initialize(params.seed);
  habitat.generate_outer(params.datapath);
  habitat.dlmwrite_outer(params.mcmcpath);
  graph.generate_grid(params.datapath,params.gridpath,
		      habitat,params.nDemes,params.nIndiv);
  graph.dlmwrite_grid(params.mcmcpath);
  o = graph.get_num_obsrv_demes();
  d = graph.get_num_total_demes();
  n = params.nIndiv;
  p = params.nSites;
  initialize_diffs();
  if (params.diploid) {
    Binvconst = 1.0; Wconst = 2.0;
  } else {
    Binvconst = 4.0; Wconst = 1.0;
  }
}
EEMS::~EEMS( ) { }
string EEMS::datapath( ) const { return params.datapath; }
string EEMS::mcmcpath( ) const { return params.mcmcpath; }
string EEMS::prevpath( ) const { return params.prevpath; }
string EEMS::gridpath( ) const { return params.gridpath; }
// Draw points randomly inside the habitat: the habitat is two-dimensional, so
// a point is represented as a row in a matrix with two columns
void EEMS::randpoint_in_habitat(MatrixXd &Seeds) {
  for (int i = 0 ; i < Seeds.rows() ; i++ ) {
    bool in = false;
    double x,y;
    while (!in) {
      x = habitat.get_xmin() + habitat.get_xspan() * draw.runif();
      y = habitat.get_ymin() + habitat.get_yspan() * draw.runif();
      in = habitat.in_point(x,y);
    }
    Seeds(i,0) = x;
    Seeds(i,1) = y;
  }
}
void EEMS::rnorm_effects(const double HalfInterval, const double rateS2, VectorXd &Effcts) {
  for (int i = 0 ; i < Effcts.rows() ; i++ ) {
    Effcts(i) = draw.rtrnorm(0.0,rateS2,HalfInterval);
  }
}
void EEMS::initialize_diffs( ) {
  cerr << "[Diffs::initialize]" << endl;
  n_2 = (double)n/2.0; nmin1 = n-1; logn = log(n);
  J = MatrixXd::Zero(n,o);
  cvec = VectorXd::Zero(o);
  for ( int i = 0 ; i < n ; i ++ ) {
    J(i,graph.get_deme_of_indiv(i)) = 1;
    cvec(graph.get_deme_of_indiv(i)) += 1;
  }
  Diffs = readMatrixXd(params.datapath + ".diffs");
  if ((Diffs.rows()!=n)||(Diffs.cols()!=n)) {
    cerr << "  Error reading dissimilarities matrix " << params.datapath + ".diffs" << endl
	 << "  Expect a " << n << "x" << n << " matrix of pairwise differences" << endl; exit(1);      
  }
  cerr << "  Loaded dissimilarities matrix from " << params.datapath + ".diffs" << endl;
  if (!isdistmat(Diffs)) {
    cerr << "  The dissimilarity matrix is not a full-rank distance matrix" << endl; exit(1);
  }
  L = MatrixXd::Constant(nmin1,n,-1.0);
  L.topRightCorner(nmin1,nmin1).setIdentity();
  JtDhatJ = MatrixXd::Zero(o,o);
  JtDobsJ = J.transpose()*Diffs*J;
  cerr << "[Diffs::initialize] Done." << endl << endl;
}
void EEMS::initialize_state( ) {
  cerr << "[EEMS::initialize_state]" << endl;
  nowdf = n;
  // Initialize the two Voronoi tessellations
  nowqtiles = 1;
  nowmtiles = draw.rnegbin(2*o,0.5); // o is the number of observed demes
  cerr << "  EEMS starts with " << nowqtiles << " qtiles and " << nowmtiles << " mtiles" << endl;
  // Draw the Voronoi centers Coord uniformly within the habitat
  nowqSeeds = MatrixXd::Zero(nowqtiles,2); randpoint_in_habitat(nowqSeeds);
  nowmSeeds = MatrixXd::Zero(nowmtiles,2); randpoint_in_habitat(nowmSeeds);
  nowmrateS2 = draw.rinvgam(0.5,0.5);
  nowqrateS2 = draw.rinvgam(0.5,0.5);
  // Assign migration rates to the Voronoi tiles
  nowqrateMu = 0;
  nowmrateMu = params.mrateMuHalfInterval*(2.0*draw.runif() - 1.0);
  // Assign rates to the Voronoi tiles
  nowqEffcts = VectorXd::Zero(nowqtiles);
  nowmEffcts = VectorXd::Zero(nowmtiles); rnorm_effects(params.mEffctHalfInterval,nowmrateS2,nowmEffcts);
  // Initialize the mapping of demes to qVoronoi tiles, i.e.,
  // For every deme in the graph -- which migration tile does the deme fall into?
  // The "color" of the deme is the index of the tile
  graph.index_closest_to_deme(nowmSeeds,nowmColors);
  // Initialize the mapping of demes to mVoronoi tiles
  graph.index_closest_to_deme(nowqSeeds,nowqColors);
  cerr << "[EEMS::initialize_state] Done." << endl << endl;
}
void EEMS::load_final_state( ) {
  cerr << "[EEMS::load_final_state]" << endl;  
  MatrixXd tempi; bool error = false;
  tempi = readMatrixXd(params.prevpath + "/lastqtiles.txt");
  if ((tempi.rows()!=1) || (tempi.cols()!=1)) { error = true; }
  nowqtiles = tempi(0,0);
  tempi = readMatrixXd(params.prevpath + "/lastmtiles.txt");
  if ((tempi.rows()!=1) || (tempi.cols()!=1)) { error = true; }
  nowmtiles = tempi(0,0);
  cerr << "  EEMS starts with " << nowqtiles << " qtiles and " << nowmtiles << " mtiles" << endl;
  tempi = readMatrixXd(params.prevpath + "/lastthetas.txt");
  if ((tempi.rows()!=1) || (tempi.cols()!=2)) { error = true; }
  nowdf = tempi(0,1);
  tempi = readMatrixXd(params.prevpath + "/lastdfpars.txt");
  if ((tempi.rows()!=1) || (tempi.cols()!=2)) { error = true; }
  params.dfmin = tempi(0,0);
  params.dfmax = tempi(0,1);
  tempi = readMatrixXd(params.prevpath + "/lastqhyper.txt");
  if ((tempi.rows()!=1) || (tempi.cols()!=1)) { error = true; }
  nowqrateS2 = tempi(0,0);
  tempi = readMatrixXd(params.prevpath + "/lastmhyper.txt");
  if ((tempi.rows()!=1) || (tempi.cols()!=2)) { error = true; }
  nowmrateMu = tempi(0,0);
  nowmrateS2 = tempi(0,1);
  tempi = readMatrixXd(params.prevpath + "/lastqeffct.txt");
  if ((tempi.rows()!=nowqtiles) || (tempi.cols()!=1)) { error = true; }
  nowqEffcts = tempi.col(0);
  tempi = readMatrixXd(params.prevpath + "/lastmeffct.txt");
  if ((tempi.rows()!=nowmtiles) || (tempi.cols()!=1)) { error = true; }
  nowmEffcts = tempi.col(0);
  nowqSeeds = readMatrixXd(params.prevpath + "/lastqseeds.txt");
  if ((nowqSeeds.rows()!=nowqtiles) || (nowqSeeds.cols()!=2)) { error = true; }
  nowmSeeds = readMatrixXd(params.prevpath + "/lastmseeds.txt");
  if ((nowmSeeds.rows()!=nowmtiles) || (nowmSeeds.cols()!=2)) { error = true; }
  // Initialize the mapping of demes to qVoronoi tiles
  graph.index_closest_to_deme(nowmSeeds,nowmColors);
  // Initialize the mapping of demes to mVoronoi tiles
  graph.index_closest_to_deme(nowqSeeds,nowqColors);
  if (error) {
    cerr << "  Error loading MCMC state from " << params.prevpath << endl; exit(1);
  }
  cerr << "[EEMS::load_final_state] Done." << endl << endl;
}
bool EEMS::start_eems(const MCMC &mcmc) {
  bool error = false;
  // The deviation of move proposals is scaled by the habitat range
  params.mSeedsProposalS2x = params.mSeedsProposalS2 * habitat.get_xspan();
  params.mSeedsProposalS2y = params.mSeedsProposalS2 * habitat.get_yspan();
  params.qSeedsProposalS2x = params.qSeedsProposalS2 * habitat.get_xspan();
  params.qSeedsProposalS2y = params.qSeedsProposalS2 * habitat.get_yspan();
  // MCMC draws are stored in memory, rather than saved to disk,
  // so it is important to thin
  int niters = mcmc.num_iters_to_save();
  mcmcmhyper = MatrixXd::Zero(niters,2);
  mcmcqhyper = MatrixXd::Zero(niters,2);
  mcmcthetas = MatrixXd::Zero(niters,2);
  mcmcpilogl = MatrixXd::Zero(niters,2);
  mcmcmtiles = VectorXd::Zero(niters);
  mcmcqtiles = VectorXd::Zero(niters);
  mcmcmRates.clear();
  mcmcqRates.clear();
  mcmcxCoord.clear();
  mcmcyCoord.clear();
  mcmcwCoord.clear();
  mcmczCoord.clear();
  nowpi = eems2_prior(nowmSeeds,nowmEffcts,nowmrateMu,nowmrateS2,
		      nowqSeeds,nowqEffcts,nowqrateMu,nowqrateS2,
		      nowdf);
  nowll = eems2_likelihood(nowmSeeds, nowmEffcts, nowmrateMu,
			   nowqSeeds, nowqEffcts, nowqrateMu,
			   nowdf);
  cerr << "Input parameters: " << endl << params << endl
       << "Initial log prior: " << nowpi << endl
       << "Initial log llike: " << nowll << endl << endl;
  if ((nowpi==-Inf) || (nowpi==Inf) || (nowll==-Inf) || (nowll==Inf)) { error = true; }
  return(error);
}
MoveType EEMS::choose_move_type( ) {
  double u1 = draw.runif( );
  double u2 = draw.runif( );
  // There are 4 types of proposals:
  // * birth/death (with equal probability)
  // * move a tile (chosen uniformly at random)
  // * update the rate of a tile (chosen uniformly at random)
  // * update the mean migration rate or the degrees of freedom (with equal probability)
  MoveType move = UNKNOWN_MOVE_TYPE;
  if (u1 < 0.25) {
    // Propose birth/death to update the Voronoi tessellation of the effective diversity,
    // with probability params.qVoronoiPr (which is 0.05 by default). Otherwise,
    // propose birth/death to update the Voronoi tessellation of the effective migration.
    move = M_VORONOI_BIRTH_DEATH;
  } else if (u1 < 0.5) {
    move = M_VORONOI_POINT_MOVE;
  } else if (u1 < 0.75) {
    move = M_VORONOI_RATE_UPDATE;
  } else {
    if (u2 < 0.5) {
      move = M_MEAN_RATE_UPDATE;
    } else {
      move = DF_UPDATE;
    }
  }
  return(move);
}
double EEMS::eval_proposal_rate_one_mtile(Proposal &proposal) const {
  return(eems2_likelihood(nowmSeeds, proposal.newmEffcts, nowmrateMu, nowqSeeds, nowqEffcts, nowqrateMu, nowdf));
}
double EEMS::eval_proposal_overall_mrate(Proposal &proposal) const {
  return(eems2_likelihood(nowmSeeds, nowmEffcts, proposal.newmrateMu, nowqSeeds, nowqEffcts, nowqrateMu, nowdf));
}
// Propose to move one tile in the migration Voronoi tessellation
double EEMS::eval_proposal_move_one_mtile(Proposal &proposal) const {
  return(eems2_likelihood(proposal.newmSeeds, nowmEffcts, nowmrateMu, nowqSeeds, nowqEffcts, nowqrateMu, nowdf));
}
double EEMS::eval_birthdeath_mVoronoi(Proposal &proposal) const {
  return(eems2_likelihood(proposal.newmSeeds, proposal.newmEffcts, nowmrateMu, nowqSeeds, nowqEffcts, nowqrateMu, nowdf));
}
void EEMS::propose_df(Proposal &proposal,const MCMC &mcmc) {
  proposal.move = DF_UPDATE;
  proposal.newpi = -Inf;
  proposal.newll = -Inf;
  // EEMS is initialized with df = nIndiv
  // Keep df = nIndiv for the first mcmc.numBurnIter/2 iterations
  // This should make it easier to move in the parameter space
  // since the likelihood is proportional to 0.5 * pdf * ll_atfixdf
  if (mcmc.currIter > (mcmc.numBurnIter/2)) {
    double newdf = draw.rnorm(nowdf,params.dfProposalS2);
    //double newdf_2 = 0.5 * newdf;
    if ( (newdf>params.dfmin) && (newdf<params.dfmax) ) {
      proposal.newdf = newdf;
      proposal.newpi = eems2_prior(nowmSeeds,nowmEffcts,nowmrateMu,nowmrateS2,
				   nowqSeeds,nowqEffcts,nowqrateMu,nowqrateS2,
				   newdf);
      proposal.newll = eems2_likelihood(nowmSeeds, nowmEffcts, nowmrateMu,
					nowqSeeds, nowqEffcts, nowqrateMu,
					newdf);
    }
  }
}
void EEMS::propose_rate_one_mtile(Proposal &proposal) {
  // Choose a tile at random to update
  int mtile = draw.runif_int(0,nowmtiles-1);
  // Make a random-walk proposal, i.e., add small offset to current value
  double curmEffct = nowmEffcts(mtile);
  double newmEffct = draw.rnorm(curmEffct,params.mEffctProposalS2);
  proposal.move = M_VORONOI_RATE_UPDATE;
  proposal.newmEffcts = nowmEffcts;
  proposal.newmEffcts(mtile) = newmEffct;
  // The prior distribution on the tile effects is truncated normal
  // So first check whether the proposed value is in range
  // Then update the prior and evaluate the new likelihood
  //  - dtrnormln(nowmEffct,0.0,nowmrateS2,params.mEffctHalfInterval) : old prior component associated with this tile
  //  + dtrnormln(newmEffct,0.0,nowmrateS2,params.mEffctHalfInterval) : new prior component associated with this tile
  if ( abs(newmEffct) < params.mEffctHalfInterval ) {
    proposal.newpi = eems2_prior(nowmSeeds,proposal.newmEffcts,nowmrateMu,nowmrateS2,
				 nowqSeeds,nowqEffcts,nowqrateMu,nowqrateS2,
				 nowdf);
    proposal.newll = eval_proposal_rate_one_mtile(proposal);
  } else {
    proposal.newpi = -Inf;
    proposal.newll = -Inf;
  }
}
void EEMS::propose_overall_mrate(Proposal &proposal) {
  // Make a random-walk Metropolis-Hastings proposal
  double newmrateMu = draw.rnorm(nowmrateMu,params.mrateMuProposalS2);
  proposal.move = M_MEAN_RATE_UPDATE;
  proposal.newmrateMu = newmrateMu;
  // If the proposed value is in range, the prior probability does not change
  // as the prior distribution on mrateMu is uniform
  // Otherwise, setting the prior and the likelihood to -Inf forces a rejection
  if ( abs(newmrateMu) < params.mrateMuHalfInterval ) {
    proposal.newpi = nowpi;
    proposal.newll = eval_proposal_overall_mrate(proposal);
  } else {
    proposal.newpi = -Inf;
    proposal.newll = -Inf;
  }
}
void EEMS::propose_move_one_mtile(Proposal &proposal) {
  // Choose a tile at random to move
  int mtile = draw.runif_int(0,nowmtiles-1);
  // Make a random-walk proposal, i.e., add small offset to current value
  // In this case, there are actually two values -- longitude and latitude
  double newmSeedx = draw.rnorm(nowmSeeds(mtile,0),params.mSeedsProposalS2x);
  double newmSeedy = draw.rnorm(nowmSeeds(mtile,1),params.mSeedsProposalS2y);
  proposal.move = M_VORONOI_POINT_MOVE;
  proposal.newmSeeds = nowmSeeds;
  proposal.newmSeeds(mtile,0) = newmSeedx;
  proposal.newmSeeds(mtile,1) = newmSeedy;
  if (habitat.in_point(newmSeedx,newmSeedy)) {
    proposal.newpi = nowpi;
    proposal.newll = eval_proposal_move_one_mtile(proposal);
  } else {
    proposal.newpi = -Inf;
    proposal.newll = -Inf;
  }
}
void EEMS::propose_birthdeath_mVoronoi(Proposal &proposal) {
  int newmtiles = nowmtiles,r;
  double u = draw.runif();
  double pBirth = 0.5;
  double pDeath = 0.5;
  proposal.newmEffcts = nowmEffcts;
  proposal.newmSeeds = nowmSeeds;
  if ((nowmtiles==1) || (u<0.5)) {  // Propose birth
    if (nowmtiles==1) { pBirth = 1.0; }
    newmtiles++;
    MatrixXd newmSeed = MatrixXd::Zero(1,2);
    randpoint_in_habitat(newmSeed);
    pairwise_distance(nowmSeeds,newmSeed).col(0).minCoeff(&r);
    double nowmEffct = nowmEffcts(r);
    double newmEffct = draw.rtrnorm(nowmEffct,params.mEffctProposalS2,params.mEffctHalfInterval);
    insertRow(proposal.newmSeeds,newmSeed.row(0));
    insertElem(proposal.newmEffcts,newmEffct);
    // Compute log(prior ratio) and log(proposal ratio)
    proposal.newratioln = log(pDeath/pBirth)
      - dtrnormln(newmEffct,nowmEffct,params.mEffctProposalS2,params.mEffctHalfInterval);

    // Two ways to get the new log prior, newpi. Should give the same answer.
    //proposal.newpi = nowpi
    //  + log((nowmtiles+params.negBiSize)/(newmtiles/params.negBiProb))
    //  + dtrnormln(newmEffct,0.0,nowmrateS2,params.mEffctHalfInterval);
    proposal.newpi = eems2_prior(proposal.newmSeeds,proposal.newmEffcts,nowmrateMu,nowmrateS2,
				 nowqSeeds,nowqEffcts,nowqrateMu,nowqrateS2,
				 nowdf);
  } else {                       // Propose death
    if (nowmtiles==2) { pBirth = 1.0; }
    newmtiles--;
    int mtileToRemove = draw.runif_int(0,newmtiles);
    MatrixXd oldmSeed = nowmSeeds.row(mtileToRemove);
    removeRow(proposal.newmSeeds,mtileToRemove);
    removeElem(proposal.newmEffcts,mtileToRemove);
    pairwise_distance(proposal.newmSeeds,oldmSeed).col(0).minCoeff(&r);
    double nowmEffct = proposal.newmEffcts(r);
    double oldmEffct = nowmEffcts(mtileToRemove);
    // Compute log(prior ratio) and log(proposal ratio)
    proposal.newratioln = log(pBirth/pDeath)
      + dtrnormln(oldmEffct,nowmEffct,params.mEffctProposalS2,params.mEffctHalfInterval);

    // Two ways to get the new log prior, newpi. Should give the same answer.
    //proposal.newpi = nowpi
    //  + log((nowmtiles/params.negBiProb)/(newmtiles+params.negBiSize))
    //  - dtrnormln(oldmEffct,0.0,nowmrateS2,params.mEffctHalfInterval);
    proposal.newpi = eems2_prior(proposal.newmSeeds,proposal.newmEffcts,nowmrateMu,nowmrateS2,
				 nowqSeeds,nowqEffcts,nowqrateMu,nowqrateS2,
				 nowdf);
  }
  proposal.move = M_VORONOI_BIRTH_DEATH;
  proposal.newmtiles = newmtiles;
  proposal.newll = eval_birthdeath_mVoronoi(proposal);
}
void EEMS::update_hyperparams( ) {
  double SSq = nowqEffcts.squaredNorm();
  double SSm = nowmEffcts.squaredNorm();
  nowqrateS2 = draw.rinvgam(params.qrateShape_2 + 0.5 * nowqtiles, params.qrateScale_2 + 0.5 * SSq);
  nowmrateS2 = draw.rinvgam(params.mrateShape_2 + 0.5 * nowmtiles, params.mrateScale_2 + 0.5 * SSm);
  nowpi = eems2_prior(nowmSeeds,nowmEffcts,nowmrateMu,nowmrateS2,
		      nowqSeeds,nowqEffcts,nowqrateMu,nowqrateS2,
		      nowdf);
}
bool EEMS::accept_proposal(Proposal &proposal) {
  double u = draw.runif( );
  // The proposal cannot be accepted because the prior is 0
  // This can happen if the proposed value falls outside the parameter's support
  if ( proposal.newpi == -Inf ) {
    proposal.newpi = nowpi;
    proposal.newll = nowll;
    return false;
  }
  double ratioln = proposal.newpi - nowpi + proposal.newll - nowll;
  // If the proposal is either birth or death, add the log(proposal ratio)
  if (proposal.move==M_VORONOI_BIRTH_DEATH) {
    ratioln += proposal.newratioln;
  }
  if ( log(u) < min(0.0,ratioln) ) {
    switch (proposal.move) {
    case M_VORONOI_RATE_UPDATE:
      nowmEffcts = proposal.newmEffcts;
      break;
    case M_MEAN_RATE_UPDATE:
      nowmrateMu = proposal.newmrateMu;
      break;
    case M_VORONOI_POINT_MOVE:
      nowmSeeds = proposal.newmSeeds;
      // Update the mapping of demes to mVoronoi tiles
      graph.index_closest_to_deme(nowmSeeds,nowmColors);
      break;
    case M_VORONOI_BIRTH_DEATH:
      nowmSeeds = proposal.newmSeeds;
      nowmEffcts = proposal.newmEffcts;
      nowmtiles = proposal.newmtiles;
      // Update the mapping of demes to mVoronoi tiles
      graph.index_closest_to_deme(nowmSeeds,nowmColors);
      break;
    case DF_UPDATE:
      nowdf = proposal.newdf;
      break;
    default:
      cerr << "[RJMCMC] Unknown move type" << endl;
      exit(1);
    }
    nowpi = proposal.newpi;
    nowll = proposal.newll;
    return true;
  } else {
    proposal.newpi = nowpi;
    proposal.newll = nowll;
    return false;
  }
}
///////////////////////////////////////////
void EEMS::print_iteration(const MCMC &mcmc) const {
  cerr << " Ending iteration " << mcmc.currIter
       << " with acceptance proportions:" << endl << mcmc
       << " and effective degrees of freedom = " << nowdf << endl
       << "         number of qVoronoi tiles = " << nowqtiles << endl
       << "         number of mVoronoi tiles = " << nowmtiles << endl
       << "          Log prior = " << nowpi << endl
       << "          Log llike = " << nowll << endl;
}
void EEMS::save_iteration(const MCMC &mcmc) {
  int iter = mcmc.to_save_iteration( );
  mcmcthetas(iter,0) = -9999;
  mcmcthetas(iter,1) = nowdf;
  mcmcqhyper(iter,0) = 0.0;
  mcmcqhyper(iter,1) = nowqrateS2;
  mcmcmhyper(iter,0) = nowmrateMu;
  mcmcmhyper(iter,1) = nowmrateS2;
  mcmcpilogl(iter,0) = nowpi;
  mcmcpilogl(iter,1) = nowll;
  mcmcqtiles(iter) = nowqtiles;
  mcmcmtiles(iter) = nowmtiles;
  for ( int t = 0 ; t < nowqtiles ; t++ ) {
    mcmcqRates.push_back(pow(10.0,nowqEffcts(t)));
  }
  for ( int t = 0 ; t < nowqtiles ; t++ ) {
    mcmcwCoord.push_back(nowqSeeds(t,0));
  }
  for ( int t = 0 ; t < nowqtiles ; t++ ) {
    mcmczCoord.push_back(nowqSeeds(t,1));
  }
  for ( int t = 0 ; t < nowmtiles ; t++ ) {
    mcmcmRates.push_back(pow(10.0,nowmEffcts(t) + nowmrateMu));
  }
  for ( int t = 0 ; t < nowmtiles ; t++ ) {
    mcmcxCoord.push_back(nowmSeeds(t,0));
  }
  for ( int t = 0 ; t < nowmtiles ; t++ ) {
    mcmcyCoord.push_back(nowmSeeds(t,1));
  }
  /////////////////////////
  MatrixXd M = MatrixXd::Zero(d,d);
  int alpha, beta;
  // Transform the log10 migration parameters into migration rates on the original scale
  for ( int edge = 0 ; edge < graph.get_num_edges() ; edge++ ) {
    graph.get_edge(edge,alpha,beta);
    double log10m_alpha = nowmEffcts(nowmColors(alpha)) + nowmrateMu;
    double log10m_beta = nowmEffcts(nowmColors(beta)) + nowmrateMu;
    M(alpha,beta) = 0.5 * pow(10.0,log10m_alpha) + 0.5 * pow(10.0,log10m_beta);
    M(beta,alpha) = M(alpha,beta);
  }  
  JtDhatJ += resistance_distance(M,o);
}
bool EEMS::output_current_state( ) const {
  ofstream out; bool error = false;
  out.open((params.mcmcpath + "/lastqtiles.txt").c_str(),ofstream::out);
  if (!out.is_open()) { error = true; return(error); }
  out << nowqtiles << endl;
  out.close( );
  out.open((params.mcmcpath + "/lastmtiles.txt").c_str(),ofstream::out);
  if (!out.is_open()) { error = true; return(error); }
  out << nowmtiles << endl;
  out.close( );
  out.open((params.mcmcpath + "/lastthetas.txt").c_str(),ofstream::out);
  if (!out.is_open()) { error = true; return(error); }
  out << fixed << setprecision(6) << -9999 << " " << nowdf << endl;
  out.close( );
  out.open((params.mcmcpath + "/lastdfpars.txt").c_str(),ofstream::out);
  if (!out.is_open()) { error = true; return(error); }
  out << fixed << setprecision(6) << params.dfmin << " " << params.dfmax << endl;
  out.close( );  
  out.open((params.mcmcpath + "/lastqhyper.txt").c_str(),ofstream::out);
  if (!out.is_open()) { error = true; return(error); }
  out << fixed << setprecision(6) << nowqrateS2 << endl;
  out.close( );
  out.open((params.mcmcpath + "/lastmhyper.txt").c_str(),ofstream::out);
  if (!out.is_open()) { error = true; return(error); }
  out << fixed << setprecision(6) << nowmrateMu << " " << nowmrateS2 << endl;
  out.close( );
  out.open((params.mcmcpath + "/lastpilogl.txt").c_str(),ofstream::out);
  if (!out.is_open()) { error = true; return(error); }
  out << fixed << setprecision(6) << nowpi << " " << nowll << endl;
  out.close( );
  out.open((params.mcmcpath + "/lastmeffct.txt").c_str(),ofstream::out);
  if (!out.is_open()) { error = true; return(error); }
  out << fixed << setprecision(6) << nowmEffcts << endl;
  out.close( );
  out.open((params.mcmcpath + "/lastmseeds.txt").c_str(),ofstream::out);
  if (!out.is_open()) { error = true; return(error); }
  out << fixed << setprecision(6) << nowmSeeds << endl;
  out.close( );
  out.open((params.mcmcpath + "/lastqeffct.txt").c_str(),ofstream::out);
  if (!out.is_open()) { error = true; return(error); }
  out << fixed << setprecision(6) << nowqEffcts << endl;
  out.close( );
  out.open((params.mcmcpath + "/lastqseeds.txt").c_str(),ofstream::out);
  if (!out.is_open()) { error = true; return(error); }
  out << fixed << setprecision(6) << nowqSeeds << endl;
  out.close( );
  return(error);
}
bool EEMS::output_results(const MCMC &mcmc) const {
  ofstream out; bool error = false;
  MatrixXd oDemes = MatrixXd::Zero(o,3);
  oDemes << graph.get_the_obsrv_demes(),cvec;
  out.open((params.mcmcpath + "/rdistoDemes.txt").c_str(),ofstream::out);
  if (!out.is_open()) { return false; }
  out << oDemes << endl;
  out.close( );
  out.open((params.mcmcpath + "/rdistJtDobsJ.txt").c_str(),ofstream::out);
  if (!out.is_open()) { return false; }
  MatrixXd Pairs = cvec*cvec.transpose(); Pairs -= cvec.asDiagonal();
  out << JtDobsJ.cwiseQuotient(Pairs) << endl;
  out.close( );
  out.open((params.mcmcpath + "/rdistJtDhatJ.txt").c_str(),ofstream::out);
  if (!out.is_open()) { return false; }
  int niters = mcmc.num_iters_to_save();
  out << JtDhatJ/niters << endl;
  out.close( );
  out.open((params.mcmcpath + "/mcmcqtiles.txt").c_str(),ofstream::out);
  if (!out.is_open()) { return false; }
  out << mcmcqtiles << endl;
  out.close( );
  out.open((params.mcmcpath + "/mcmcmtiles.txt").c_str(),ofstream::out);
  if (!out.is_open()) { return false; }
  out << mcmcmtiles << endl;
  out.close( );
  out.open((params.mcmcpath + "/mcmcthetas.txt").c_str(),ofstream::out);
  if (!out.is_open()) { return false; }
  out << fixed << setprecision(6) << mcmcthetas << endl;
  out.close( );
  out.open((params.mcmcpath + "/mcmcqhyper.txt").c_str(),ofstream::out);
  if (!out.is_open()) { return false; }
  out << fixed << setprecision(6) << mcmcqhyper << endl;
  out.close( );
  out.open((params.mcmcpath + "/mcmcmhyper.txt").c_str(),ofstream::out);
  if (!out.is_open()) { return false; }
  out << fixed << setprecision(6) << mcmcmhyper << endl;
  out.close( );
  out.open((params.mcmcpath + "/mcmcpilogl.txt").c_str(),ofstream::out);
  if (!out.is_open()) { return false; }
  out << fixed << setprecision(6) << mcmcpilogl << endl;
  out.close( );
  error = dlmcell(params.mcmcpath + "/mcmcmrates.txt",mcmcmtiles,mcmcmRates); if (error) { return(error); }
  error = dlmcell(params.mcmcpath + "/mcmcxcoord.txt",mcmcmtiles,mcmcxCoord); if (error) { return(error); }
  error = dlmcell(params.mcmcpath + "/mcmcycoord.txt",mcmcmtiles,mcmcyCoord); if (error) { return(error); }
  error = dlmcell(params.mcmcpath + "/mcmcqrates.txt",mcmcqtiles,mcmcqRates); if (error) { return(error); }
  error = dlmcell(params.mcmcpath + "/mcmcwcoord.txt",mcmcqtiles,mcmcwCoord); if (error) { return(error); }
  error = dlmcell(params.mcmcpath + "/mcmczcoord.txt",mcmcqtiles,mcmczCoord); if (error) { return(error); }
  error = output_current_state( ); if (error) { return(error); }
  out.open((params.mcmcpath + "/eemsrun.txt").c_str(),ofstream::out);
  if (!out.is_open( )) { return false; }
  out << "Input parameter values:" << endl << params << endl
      << "Acceptance proportions:" << endl << mcmc << endl
      << "Final log prior: " << nowpi << endl
      << "Final log llike: " << nowll << endl;
  out.close( );
  cerr << "Final log prior: " << nowpi << endl
       << "Final log llike: " << nowll << endl;
  return(error);
}
double EEMS::eems2_prior(const MatrixXd &mSeeds, const VectorXd &mEffcts, const double mrateMu, const double mrateS2,
			 const MatrixXd &qSeeds, const VectorXd &qEffcts, const double qrateMu, const double qrateS2,
			 const double df) const {
  bool inrange = true;
  int qtiles = qEffcts.size();
  int mtiles = mEffcts.size();
  // First check that all parameters fall into their support range
  for ( int i = 0 ; i < qtiles ; i++ ) {
    if (!habitat.in_point(qSeeds(i,0),qSeeds(i,1))) { inrange = false; }
  }
  for ( int i = 0 ; i < mtiles ; i++ ) {
    if (!habitat.in_point(mSeeds(i,0),mSeeds(i,1))) { inrange = false; }
  }
  if (qEffcts.cwiseAbs().minCoeff()>params.qEffctHalfInterval) { inrange = false; }
  if (mEffcts.cwiseAbs().minCoeff()>params.mEffctHalfInterval) { inrange = false; }
  if (abs(mrateMu)>params.mrateMuHalfInterval) { inrange = false; }
  if ((df<params.dfmin) || (df>params.dfmax)) { inrange = false; }
  if (!inrange) { return (-Inf); }
  // Then compute the prior, on the log scale
  double logpi = - log(df)
    + dnegbinln(mtiles,params.negBiSize,params.negBiProb)
    + dnegbinln(qtiles,params.negBiSize,params.negBiProb)
    + dinvgamln(mrateS2,params.mrateShape_2,params.mrateScale_2)
    + dinvgamln(qrateS2,params.qrateShape_2,params.qrateScale_2);
  for (int i = 0 ; i < qtiles ; i++) {
    logpi += dtrnormln(qEffcts(i),0.0,qrateS2,params.qEffctHalfInterval);
  }
  for (int i = 0 ; i < mtiles ; i++) {
    logpi += dtrnormln(mEffcts(i),0.0,mrateS2,params.mEffctHalfInterval);
  }
  return (logpi);
}
double EEMS::eems2_likelihood(const MatrixXd &mSeeds, const VectorXd &mEffcts, const double mrateMu,
			      const MatrixXd &qSeeds, const VectorXd &qEffcts, const double qrateMu,
			      const double df) const {

  double sigma2 = 1.0;
  
  // mSeeds, mEffcts and mrateMu define the migration Voronoi tessellation
  // qSeeds, qEffcts define the diversity Voronoi tessellation
  // These are EEMS parameters, so no need to pass them to test_likelihood
  VectorXi mColors, qColors;
  // For every deme in the graph -- which migration tile does the deme fall into? 
  graph.index_closest_to_deme(mSeeds,mColors);
  // For every deme in the graph -- which diversity tile does the deme fall into? 
  graph.index_closest_to_deme(qSeeds,qColors);
  VectorXd W = VectorXd::Zero(d);
  // Transform the log10 diversity parameters into diversity rates on the original scale
  for ( int alpha = 0 ; alpha < d ; alpha++ ) {
    double log10q_alpha = qEffcts(qColors(alpha)); // qrateMu = 0.0
    W(alpha) = pow(10.0,log10q_alpha);
  }
  MatrixXd M = MatrixXd::Zero(d,d);
  int alpha, beta;
  // Transform the log10 migration parameters into migration rates on the original scale
  for ( int edge = 0 ; edge < graph.get_num_edges() ; edge++ ) {
    graph.get_edge(edge,alpha,beta);
    double log10m_alpha = mEffcts(mColors(alpha)) + mrateMu;
    double log10m_beta = mEffcts(mColors(beta)) + mrateMu;
    M(alpha,beta) = 0.5 * pow(10.0,log10m_alpha) + 0.5 * pow(10.0,log10m_beta);
    M(beta,alpha) = M(alpha,beta);
  }
  // J is an indicator matrix such that J(i,a) = 1 if individual i comes from deme a,
  // and J(i,a) = 0 otherwise  
  MatrixXd Delta = expected_between_dissimilarities(J, M * Binvconst);
  double logll = pseudowishpdfln(-L*Diffs*L.transpose(),
				 -L*Delta*L.transpose()*sigma2/df,df);
  return (logll);
}
