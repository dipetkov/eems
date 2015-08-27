
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
  cinv = VectorXd::Zero(o);
  for ( int i = 0 ; i < n ; i ++ ) {
    J(i,graph.get_deme_of_indiv(i)) = 1;
    cvec(graph.get_deme_of_indiv(i)) += 1;
  }
  cinv = pow(cvec.array(),-1.0).matrix();  // cinv is the vector of inverse counts
  cmin1 = cvec; cmin1.array() -= 1;        // cmin1 is the vector of counts - 1
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
  ldLLt = logdet(L*L.transpose());
  ldLDLt = logdet(-L*Diffs*L.transpose());
  ldDiQ = ldLLt - ldLDLt;
  cerr << "[Diffs::initialize] Done." << endl << endl;
}
void EEMS::initialize_state( ) {
  cerr << "[EEMS::initialize_state]" << endl;
  nowdf = n;
  // Initialize the two Voronoi tessellations
  nowsigma2 = draw.rinvgam(3.0,1.0);
  nowqtiles = draw.rnegbin(2*o,0.5); // o is the number of observed demes
  nowmtiles = draw.rnegbin(2*o,0.5);
  cerr << "  EEMS starts with " << nowqtiles << " qtiles and " << nowmtiles << " mtiles" << endl;
  // Draw the Voronoi centers Coord uniformly within the habitat
  nowqSeeds = MatrixXd::Zero(nowqtiles,2); randpoint_in_habitat(nowqSeeds);
  nowmSeeds = MatrixXd::Zero(nowmtiles,2); randpoint_in_habitat(nowmSeeds);
  nowmrateS2 = draw.rinvgam(0.5,0.5);
  nowqrateS2 = draw.rinvgam(0.5,0.5);
  // Assign migration rates to the Voronoi tiles
  nowmrateMu = params.mrateMuHalfInterval*(2.0*draw.runif() - 1.0);
  // Assign rates to the Voronoi tiles
  nowqEffcts = VectorXd::Zero(nowqtiles); rnorm_effects(params.qEffctHalfInterval,nowqrateS2,nowqEffcts);
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
  nowsigma2 = tempi(0,0);
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
  this->eval_prior();
  this->eval_likelihood();
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
    if (u2 < params.qVoronoiPr) {
      move = Q_VORONOI_BIRTH_DEATH;
    } else {
      move = M_VORONOI_BIRTH_DEATH;
    }
  } else if (u1 < 0.5) {
    if (u2 < params.qVoronoiPr) {
      move = Q_VORONOI_POINT_MOVE;
    } else {
      move = M_VORONOI_POINT_MOVE;
    }
  } else if (u1 < 0.75) {
    if (u2 < params.qVoronoiPr) {
      move = Q_VORONOI_RATE_UPDATE;
    } else {
      move = M_VORONOI_RATE_UPDATE;
    }
  } else {
    if (u2 < 0.5) {
      move = M_MEAN_RATE_UPDATE;
    } else {
      move = DF_UPDATE;
    }
  }
  return(move);
}
double EEMS::eval_prior( ) {
  // The parameters should always be in range
  bool inrange = true;
  for ( int i = 0 ; i < nowqtiles ; i++ ) {
    if (!habitat.in_point(nowqSeeds(i,0),nowqSeeds(i,1))) { inrange = false; }
  }
  for ( int i = 0 ; i < nowmtiles ; i++ ) {
    if (!habitat.in_point(nowmSeeds(i,0),nowmSeeds(i,1))) { inrange = false; }
  }
  if (nowqEffcts.cwiseAbs().minCoeff()>params.qEffctHalfInterval) { inrange = false; }
  if (nowmEffcts.cwiseAbs().minCoeff()>params.mEffctHalfInterval) { inrange = false; }
  if (abs(nowmrateMu)>params.mrateMuHalfInterval) { inrange = false; }
  if ((nowdf<params.dfmin) || (nowdf>params.dfmax)) { inrange = false; }
  if (!inrange) { return (-Inf); }
  // Use the normal pdf with mean mu and standard deviation sigma to compute
  // normalizing constant for a truncated normal pdf with mean mu, standard
  // deviation sigma and support [lower bound, upper bound].
  boost::math::normal nowmrateNorm(0.0,sqrt(nowmrateS2));
  boost::math::normal nowqrateNorm(0.0,sqrt(nowqrateS2));
  nowpi = - log(nowdf)
    + lgamma(params.negBiSize+nowmtiles) - lgamma(nowmtiles+1.0) + nowmtiles*log(params.negBiProb)
    + lgamma(params.negBiSize+nowqtiles) - lgamma(nowqtiles+1.0) + nowqtiles*log(params.negBiProb)
    - (params.mrateShape_2+1.0)*log(nowmrateS2) - params.mrateScale_2/nowmrateS2
    - (params.qrateShape_2+1.0)*log(nowqrateS2) - params.qrateScale_2/nowqrateS2
    - (params.sigmaShape_2+1.0)*log(nowsigma2) - params.sigmaScale_2/nowsigma2
    - (nowmtiles/2.0)*log(nowmrateS2) - nowmEffcts.squaredNorm()/(2.0*nowmrateS2)
    - (nowqtiles/2.0)*log(nowqrateS2) - nowqEffcts.squaredNorm()/(2.0*nowqrateS2)
    - nowmtiles * log(cdf(nowmrateNorm,params.mEffctHalfInterval) - cdf(nowmrateNorm,-params.mEffctHalfInterval))
    - nowqtiles * log(cdf(nowqrateNorm,params.qEffctHalfInterval) - cdf(nowqrateNorm,-params.qEffctHalfInterval));
  return (nowpi);
}
double EEMS::eval_likelihood( ) {
  // Expected genetic dissimilarities Delta are modeled as
  // Delta(a,b) = BetweenDistance(a,b) + ( WithinDiversity(a) + WithinDiversity(b) )/2
  //            = B(a,b) + ( W(a) + W(b) )/ 2
  // For every deme in the graph -- what is its effective diversity q(a)?
  calc_within(nowqColors,nowqEffcts,nowW);
  // For every pair of demes -- what is the effective resistance distance B(a,b)?
  // Binv is the inverse of B
  calc_between(nowmColors,nowmEffcts,nowmrateMu,nowBinv);
  double triDeltaQD, ll_atfixdf;
  // Compute the Wishart log likelihood
  nowll = EEMS_wishpdfln(nowBinv,nowW,nowsigma2,nowdf,nowtriDeltaQD,nowll_atfixdf);
  return (nowll);
}
void EEMS::calc_within(const VectorXi &qColors, const VectorXd &qEffcts, VectorXd &W) const {
  // o is the number of observed demes in the graph
  if (W.size()!=o) { W.resize(o); }
  // For every observed deme in the graph
  for ( int alpha = 0 ; alpha < o ; alpha++ ) {
    // WithinDiversity = W(a) = 10^( q_alpha ) = 10^( q_tile(tile_alpha))
    W(alpha) = pow(10.0,qEffcts(qColors(alpha)));
  }
  // The constant is slightly different for haploid and diploid species
  W *= Wconst;
}
void EEMS::calc_between(const VectorXi &mColors, const VectorXd &mEffcts, const double mrateMu, MatrixXd &Binv) const {
  // d is the number of demes in the graph (observed or not)
  if (Binv.rows()!=d||Binv.cols()!=d) { Binv.resize(d,d); }
  // A sparse matrix of migration rates will be constructed on the fly
  // First a triplet (a,b,m) is added for each edge (a,b)
  // I have decided to construct the sparse matrix rather than update it
  // because a single change to the migration Voronoi tessellation can
  // change the migration rate of many edges simultaneously
  vector<Tri> coefficients;
  int alpha, beta;
  // For every edge in the graph -- it does not matter whether demes are observed or not
  // as the resistance distance takes into consideration all paths between a and b
  // to produces the effective resistance B(a,b)
  for ( int edge = 0 ; edge < graph.get_num_edges() ; edge++ ) {
    graph.get_edge(edge,alpha,beta);
    // On the log10 scale, log10(m_alpha) = mrateMu + m_alpha = mrateMu + m_tile(tile_alpha)
    // On the log10 scale, log10(m_beta) = mrateMu + m_beta = mrateMu + m_tile(tile_beta)
    double log10m_alpha = mrateMu + mEffcts(mColors(alpha));
    double log10m_beta = mrateMu + mEffcts(mColors(beta));
    // Then on the original scale, m(alpha,beta) = (10^m_alpha + 10^m_beta)/2
    double m_ab = 0.5 * pow(10.0,log10m_alpha) + 0.5 * pow(10.0,log10m_beta);
    // The graph is undirected, so m(alpha->beta) = m(beta->alpha)
    coefficients.push_back(Tri(alpha,beta,m_ab));
    coefficients.push_back(Tri(beta,alpha,m_ab));
  }
  SpMat sparseM(d,d);
  // Actually construct and fill in the sparse matrix
  sparseM.setFromTriplets(coefficients.begin(),coefficients.end());
  // Initialize a dense matrix from the sparse matrix
  MatrixXd M = MatrixXd(sparseM);
  // Read S1.4 Computing the resistance distances in the Supplementary
  // This computation is specific to the resistance distance metric
  // but the point is that we have a method to compute the between
  // demes component of the expected genetic dissimilarities from
  // the sparse matrix of migration rates
  // Here instead of B, we compute its inverse Binv
  MatrixXd Hinv = - M; Hinv.diagonal() += M.rowwise().sum(); Hinv.array() += 1.0;
  if (o==d) {
    Binv = -0.5 * Hinv; // If all the demes are observed, it is quite easy to compute Binv
  } else {
    Binv = -0.5 * Hinv.topLeftCorner(o,o);
    Binv += 0.5 * Hinv.topRightCorner(o,d-o) *
      Hinv.bottomRightCorner(d-o,d-o).selfadjointView<Lower>().llt().solve(Hinv.bottomLeftCorner(d-o,o));
  }
  // The constant is slightly different for haploid and diploid species 
  Binv *= Binvconst;
}
/*
  This function implements the computations described in the Section S1.3 in the Supplementary Information,
  "Computing the Wishart log likelihood l(k, m, q, sigma2)", and I have tried to used similar notation
  For example, MatrixXd X = lu.solve(T) is equation (S20)
  since lu is the decomposition of (B*C - W) and T is B * inv(W)
  Returns wishpdfln( -L*D*L' ; - (sigma2/df) * L*Delta(m,q)*L' , df )
 */
double EEMS::EEMS_wishpdfln(const MatrixXd &Binv, const VectorXd &W, const double sigma2, const double df,
			    double &triDeltaQD, double &ll_atfixdf) const {
  double df_2 = 0.5 * df;
  MatrixXd T = Binv.selfadjointView<Lower>().ldlt().solve(MatrixXd::Identity(o,o));
  T *= cvec.asDiagonal(); T -= W.asDiagonal();             // Now T = B*C - W
  PartialPivLU<MatrixXd> lu(T);
  VectorXd Winv = pow(W.array(),-1.0);
  T.noalias() = T*Winv.asDiagonal()*cinv.asDiagonal();
  T += cinv.asDiagonal();                                  // Now T = B*inv(W)
  MatrixXd X = lu.solve(T);
  VectorXd Xc_Winv = X*cvec - Winv;
  double oDinvo = cvec.dot(Xc_Winv);                       // oDinvo = 1'*inv(Delta)*1
  double oDiDDi = Xc_Winv.transpose()*JtDobsJ*Xc_Winv;     // oDiDDi = 1'*inv(Delta)*D*inv(D)*1
  triDeltaQD = trace_AxB(X,JtDobsJ) - oDiDDi/oDinvo;       // triDeltaQD = tr(inv(Delta)*Q*D)
  double ldetDinvQ = logn - log(abs(oDinvo))               // ldetDinvQ = logDet(-inv(Delta)*Q)
    + cmin1.dot(Winv.array().log().matrix())               // log(det(W_n) * det(inv(W_o)))
    - lu.matrixLU().diagonal().array().abs().log().sum();  // log(abs(det(B*C - W)))
  ll_atfixdf = ldetDinvQ - triDeltaQD/sigma2 - nmin1*log(sigma2) - ldDiQ;
  return (df_2*ll_atfixdf + nmin1*df_2*log(df_2) - mvgammaln(df_2,nmin1) - n_2*ldLDLt);
}
double EEMS::eval_proposal_rate_one_qtile(Proposal &proposal) const {
  calc_within(nowqColors,proposal.newqEffcts,proposal.newW);
  return (EEMS_wishpdfln(nowBinv,proposal.newW,nowsigma2,nowdf,proposal.newtriDeltaQD,proposal.newll_atfixdf));
}
double EEMS::eval_proposal_move_one_qtile(Proposal &proposal) const {
  graph.index_closest_to_deme(proposal.newqSeeds,proposal.newqColors);
  calc_within(proposal.newqColors,nowqEffcts,proposal.newW);
  return (EEMS_wishpdfln(nowBinv,proposal.newW,nowsigma2,nowdf,proposal.newtriDeltaQD,proposal.newll_atfixdf));
}
double EEMS::eval_birthdeath_qVoronoi(Proposal &proposal) const {
  graph.index_closest_to_deme(proposal.newqSeeds,proposal.newqColors);
  calc_within(proposal.newqColors,proposal.newqEffcts,proposal.newW);
  return (EEMS_wishpdfln(nowBinv,proposal.newW,nowsigma2,nowdf,proposal.newtriDeltaQD,proposal.newll_atfixdf));
}
double EEMS::eval_proposal_rate_one_mtile(Proposal &proposal) const {
  calc_between(nowmColors,proposal.newmEffcts,nowmrateMu,proposal.newBinv);
  return (EEMS_wishpdfln(proposal.newBinv,nowW,nowsigma2,nowdf,proposal.newtriDeltaQD,proposal.newll_atfixdf));
}
double EEMS::eval_proposal_overall_mrate(Proposal &proposal) const {
  calc_between(nowmColors,nowmEffcts,proposal.newmrateMu,proposal.newBinv);
  return (EEMS_wishpdfln(proposal.newBinv,nowW,nowsigma2,nowdf,proposal.newtriDeltaQD,proposal.newll_atfixdf));
}
// Propose to move one tile in the migration Voronoi tessellation
double EEMS::eval_proposal_move_one_mtile(Proposal &proposal) const {
  graph.index_closest_to_deme(proposal.newmSeeds,proposal.newmColors);
  calc_between(proposal.newmColors,nowmEffcts,nowmrateMu,proposal.newBinv);
  return (EEMS_wishpdfln(proposal.newBinv,nowW,nowsigma2,nowdf,proposal.newtriDeltaQD,proposal.newll_atfixdf));
}
double EEMS::eval_birthdeath_mVoronoi(Proposal &proposal) const {
  graph.index_closest_to_deme(proposal.newmSeeds,proposal.newmColors);
  calc_between(proposal.newmColors,proposal.newmEffcts,nowmrateMu,proposal.newBinv);
  return (EEMS_wishpdfln(proposal.newBinv,nowW,nowsigma2,nowdf,proposal.newtriDeltaQD,proposal.newll_atfixdf));
}
///////////////////////////////////////////
void EEMS::update_sigma2( ) {
  double nowdf_2 = 0.5 * nowdf;
  nowll_atfixdf += nowtriDeltaQD/nowsigma2 + nmin1*log(nowsigma2);
  nowpi += (params.sigmaShape_2+1.0)*log(nowsigma2) + params.sigmaScale_2/nowsigma2;
  nowsigma2 = draw.rinvgam( params.sigmaShape_2 + nowdf_2*nmin1,
			    params.sigmaScale_2 + nowdf_2*nowtriDeltaQD );
  nowll_atfixdf -= nowtriDeltaQD/nowsigma2 + nmin1*log(nowsigma2);
  nowpi -= (params.sigmaShape_2+1.0)*log(nowsigma2) + params.sigmaScale_2/nowsigma2;
  nowll = nowdf_2 * nowll_atfixdf + nmin1*nowdf_2*log(nowdf_2) - mvgammaln(nowdf_2,nmin1) - n_2*ldLDLt;
}
void EEMS::propose_df(Proposal &proposal,const MCMC &mcmc) {
  proposal.move = DF_UPDATE;
  proposal.newtriDeltaQD = nowtriDeltaQD;
  proposal.newll_atfixdf = nowll_atfixdf;
  proposal.newpi = -Inf;
  proposal.newll = -Inf;
  // EEMS is initialized with df = nIndiv
  // Keep df = nIndiv for the first mcmc.numBurnIter/2 iterations
  // This should make it easier to move in the parameter space
  // since the likelihood is proportional to 0.5 * pdf * ll_atfixdf
  if (mcmc.currIter > (mcmc.numBurnIter/2)) {
    double newdf = draw.rnorm(nowdf,params.dfProposalS2);
    double newdf_2 = 0.5 * newdf;
    if ( (newdf>params.dfmin) && (newdf<params.dfmax) ) {
      proposal.newdf = newdf;
      proposal.newpi = nowpi + log(nowdf) - log(newdf);
      proposal.newll =
	newdf_2*nowll_atfixdf + nmin1*newdf_2*log(newdf_2) - mvgammaln(newdf_2,nmin1) - n_2*ldLDLt;
    }
  }
}
void EEMS::propose_rate_one_qtile(Proposal &proposal) {
  // Choose a tile at random to update
  int qtile = draw.runif_int(0,nowqtiles-1);
  // Make a random-walk proposal, i.e., add small offset to current value
  double nowqEffct = nowqEffcts(qtile);
  double newqEffct = draw.rnorm(nowqEffct,params.qEffctProposalS2);
  proposal.move = Q_VORONOI_RATE_UPDATE;
  proposal.newqEffcts = nowqEffcts;
  proposal.newqEffcts(qtile) = newqEffct;
  // The prior distribution on the tile effects is truncated normal
  // So first check whether the proposed value is in range
  // Then update the prior and evaluate the new likelihood
  //  - dtrnormln(nowqEffct,0.0,nowqrateS2,params.qEffctHalfInterval) : old prior component associated with this tile
  //  + dtrnormln(newqEffct,0.0,nowqrateS2,params.qEffctHalfInterval) : new prior component associated with this tile
  if ( abs(newqEffct) < params.qEffctHalfInterval ) {
    proposal.newpi = nowpi
      - dtrnormln(nowqEffct,0.0,nowqrateS2,params.qEffctHalfInterval)
      + dtrnormln(newqEffct,0.0,nowqrateS2,params.qEffctHalfInterval);
    proposal.newll = eval_proposal_rate_one_qtile(proposal);
  } else {
    proposal.newpi = -Inf;
    proposal.newll = -Inf;
  }
}
void EEMS::propose_rate_one_mtile(Proposal &proposal) {
  // Choose a tile at random to update
  int mtile = draw.runif_int(0,nowmtiles-1);
  // Make a random-walk proposal, i.e., add small offset to current value
  double nowmEffct = nowmEffcts(mtile);
  double newmEffct = draw.rnorm(nowmEffct,params.mEffctProposalS2);
  proposal.move = M_VORONOI_RATE_UPDATE;
  proposal.newmEffcts = nowmEffcts;
  proposal.newmEffcts(mtile) = newmEffct;
  if ( abs(newmEffct) < params.mEffctHalfInterval ) {
    proposal.newpi = nowpi
      - dtrnormln(nowmEffct,0.0,nowmrateS2,params.mEffctHalfInterval)
      + dtrnormln(newmEffct,0.0,nowmrateS2,params.mEffctHalfInterval);
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
void EEMS::propose_move_one_qtile(Proposal &proposal) {
  // Choose a tile at random to move
  int qtile = draw.runif_int(0,nowqtiles-1);
  // Make a random-walk proposal, i.e., add small offset to current value
  // In this case, there are actually two values -- longitude and latitude
  double newqSeedx = draw.rnorm(nowqSeeds(qtile,0),params.qSeedsProposalS2x);
  double newqSeedy = draw.rnorm(nowqSeeds(qtile,1),params.qSeedsProposalS2y);
  proposal.move = Q_VORONOI_POINT_MOVE;
  proposal.newqSeeds = nowqSeeds;
  proposal.newqSeeds(qtile,0) = newqSeedx;
  proposal.newqSeeds(qtile,1) = newqSeedy;
  if (habitat.in_point(newqSeedx,newqSeedy)) {
    proposal.newpi = nowpi;
    proposal.newll = eval_proposal_move_one_qtile(proposal);
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
void EEMS::propose_birthdeath_qVoronoi(Proposal &proposal) {
  int newqtiles = nowqtiles,r;
  double u = draw.runif();
  double pBirth = 0.5;
  double pDeath = 0.5;
  proposal.newqEffcts = nowqEffcts;
  proposal.newqSeeds = nowqSeeds;
  // If there is exactly one tile, rule out a death proposal
  if ((nowqtiles==1) || (u<0.5)) { // Propose birth
    if (nowqtiles==1) { pBirth = 1.0; }
    newqtiles++;
    MatrixXd newqSeed = MatrixXd::Zero(1,2);
    randpoint_in_habitat(newqSeed);
    pairwise_distance(nowqSeeds,newqSeed).col(0).minCoeff(&r);
    // The new tile is assigned a rate by perturbing the current rate at the new seed    
    double nowqEffct = nowqEffcts(r);
    double newqEffct = draw.rtrnorm(nowqEffct,params.qEffctProposalS2,params.qEffctHalfInterval);
    insertRow(proposal.newqSeeds,newqSeed.row(0));
    insertElem(proposal.newqEffcts,newqEffct);
    // Compute log(proposal ratio) and log(prior ratio)
    proposal.newratioln = log(pDeath/pBirth)
      - dtrnormln(newqEffct,nowqEffct,params.qEffctProposalS2,params.qEffctHalfInterval);
    proposal.newpi = nowpi
      + log((nowqtiles+params.negBiSize)/(newqtiles/params.negBiProb))
      + dtrnormln(newqEffct,0.0,nowqrateS2,params.qEffctHalfInterval);
  } else {                      // Propose death
    if (nowqtiles==2) { pBirth = 1.0; }
    newqtiles--;
    int qtileToRemove = draw.runif_int(0,newqtiles);
    MatrixXd oldqSeed = nowqSeeds.row(qtileToRemove);
    removeRow(proposal.newqSeeds,qtileToRemove);
    removeElem(proposal.newqEffcts,qtileToRemove);
    pairwise_distance(proposal.newqSeeds,oldqSeed).col(0).minCoeff(&r);
    double nowqEffct = proposal.newqEffcts(r);
    double oldqEffct = nowqEffcts(qtileToRemove);
    // Compute log(prior ratio) and log(proposal ratio)
    proposal.newratioln = log(pBirth/pDeath)
      + dtrnormln(oldqEffct,nowqEffct,params.qEffctProposalS2,params.qEffctHalfInterval);
    proposal.newpi = nowpi
      + log((nowqtiles/params.negBiProb)/(newqtiles+params.negBiSize))
      - dtrnormln(oldqEffct,0.0,nowqrateS2,params.qEffctHalfInterval);
  }
  proposal.move = Q_VORONOI_BIRTH_DEATH;
  proposal.newqtiles = newqtiles;
  proposal.newll = eval_birthdeath_qVoronoi(proposal);
}
void EEMS::propose_birthdeath_mVoronoi(Proposal &proposal) {
  int newmtiles = nowmtiles,r;
  double u = draw.runif();
  double pBirth = 0.5;
  double pDeath = 0.5;
  //boost::math::normal pnorm(0.0,sqrt(nowmrateS2));
  proposal.newmEffcts = nowmEffcts;
  proposal.newmSeeds = nowmSeeds;
  if ((nowmtiles==1) || (u<0.5)) { // Propose birth
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
    proposal.newpi = nowpi
      + log((nowmtiles+params.negBiSize)/(newmtiles/params.negBiProb))
      + dtrnormln(newmEffct,0.0,nowmrateS2,params.mEffctHalfInterval);
  } else {                      // Propose death
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
    proposal.newpi = nowpi
      + log((nowmtiles/params.negBiProb)/(newmtiles+params.negBiSize))
      - dtrnormln(oldmEffct,0.0,nowmrateS2,params.mEffctHalfInterval);
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
  this->eval_prior();
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
  if (proposal.move==Q_VORONOI_BIRTH_DEATH || proposal.move==M_VORONOI_BIRTH_DEATH) {
    ratioln += proposal.newratioln;
  }
  if ( log(u) < min(0.0,ratioln) ) {
    switch (proposal.move) {
    case Q_VORONOI_RATE_UPDATE:
      nowqEffcts = proposal.newqEffcts;
      nowW = proposal.newW;
      break;
    case Q_VORONOI_POINT_MOVE:
      nowqSeeds = proposal.newqSeeds;
      nowqColors = proposal.newqColors;
      nowW = proposal.newW;
      break;
    case Q_VORONOI_BIRTH_DEATH:
      nowqSeeds = proposal.newqSeeds;
      nowqEffcts = proposal.newqEffcts;
      nowqtiles = proposal.newqtiles;
      nowqColors = proposal.newqColors;
      nowW = proposal.newW;
      break;
    case M_VORONOI_RATE_UPDATE:
      nowmEffcts = proposal.newmEffcts;
      nowBinv = proposal.newBinv;
      break;
    case M_MEAN_RATE_UPDATE:
      nowmrateMu = proposal.newmrateMu;
      nowBinv = proposal.newBinv;
      break;
    case M_VORONOI_POINT_MOVE:
      nowmSeeds = proposal.newmSeeds;
      nowmColors = proposal.newmColors;
      nowBinv = proposal.newBinv;
      break;
    case M_VORONOI_BIRTH_DEATH:
      nowmSeeds = proposal.newmSeeds;
      nowmEffcts = proposal.newmEffcts;
      nowmtiles = proposal.newmtiles;
      nowmColors = proposal.newmColors;
      nowBinv = proposal.newBinv;
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
    nowtriDeltaQD = proposal.newtriDeltaQD;
    nowll_atfixdf = proposal.newll_atfixdf;
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
  mcmcthetas(iter,0) = nowsigma2;
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
  MatrixXd B = nowBinv.inverse();
  VectorXd h = B.diagonal();    // If B = -2H, then diag(B) = -2diag(H) = -2h
  B -= 0.5 * h.replicate(1,o);  // Therefore 1h' + h1' - 2H = -1diag(B)'/2 - diag(B)1'/2 + B
  B -= 0.5 * h.transpose().replicate(o,1);
  B += 0.5 * nowW.replicate(1,o);
  B += 0.5 * nowW.transpose().replicate(o,1);
  JtDhatJ += nowsigma2 * B;
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
  out << fixed << setprecision(6) << nowsigma2 << " " << nowdf << endl;
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
void EEMS::check_ll_computation( ) const {
  double pi0 = test_prior(nowmSeeds,nowmEffcts,nowmrateMu,nowqSeeds,nowqEffcts,nowdf,nowsigma2,nowmrateS2,nowqrateS2);
  double ll0 = test_likelihood(nowmSeeds,nowmEffcts,nowmrateMu,nowqSeeds,nowqEffcts,nowdf,nowsigma2);
  if (!isfinite(nowpi) || !isfinite(pi0) ||
      !isfinite(nowll) || !isfinite(ll0) ||
      (abs(nowpi-pi0)/abs(pi0)>1e-12) ||
      (abs(nowll-ll0)/abs(ll0)>1e-12)) {
    cerr << "[EEMS::testing]   |ll0-ll|/|ll0| = " << abs(nowll - ll0)/abs(ll0) << endl;
    cerr << "[EEMS::testing]   |pi0-pi|/|pi0| = " << abs(nowpi - pi0)/abs(pi0) << endl;
    exit(1);
  }
}
double EEMS::test_prior(const MatrixXd &mSeeds, const VectorXd &mEffcts, const double mrateMu,
			const MatrixXd &qSeeds, const VectorXd &qEffcts,
			const double df, const double sigma2, const double mrateS2, const double qrateS2) const {
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
    + dinvgamln(qrateS2,params.qrateShape_2,params.qrateScale_2)
    + dinvgamln( sigma2,params.sigmaShape_2,params.sigmaScale_2);
  for (int i = 0 ; i < qtiles ; i++) {
    logpi += dtrnormln(qEffcts(i),0.0,qrateS2,params.qEffctHalfInterval);
  }
  for (int i = 0 ; i < mtiles ; i++) {
    logpi += dtrnormln(mEffcts(i),0.0,mrateS2,params.mEffctHalfInterval);
  }
  return (logpi);
}
double EEMS::test_likelihood(const MatrixXd &mSeeds, const VectorXd &mEffcts, const double mrateMu,
			     const MatrixXd &qSeeds, const VectorXd &qEffcts,
			     const double df, const double sigma2) const {
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
  MatrixXd Delta = expected_dissimilarities(J, M * Binvconst, W * Wconst);
  // Exactly equation S13
  double logll = wishpdfln(-L*Diffs*L.transpose(),
			   -L*Delta*L.transpose()*sigma2/df,df);
  return (logll);
}
