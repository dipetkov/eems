
#include "eems.hpp"

EEMS::EEMS(const Params &params,const long seed) {
  this->params = params;
  draw.initialize(seed);
  cerr << "[runEEMS] Using Boost version: ";
  get_boost_version(cerr);
  cerr << ", and Eigen version: ";
  get_eigen_version(cerr);
  cerr << endl << "          EEMS was tested with Boost 1_57 and Eigen 3.2.2" << endl << endl;
  if (!numeric_limits<double>::has_infinity) {
    cerr << "[runEEMS] numeric_limits<double>::has_infinity : no" << endl;
    exit(1);
  }
  boost::filesystem::path dir(params.mcmcpath.c_str());
  if (!boost::filesystem::create_directory(dir)) {
    cerr << "[runEEMS] Failed to create output directory " << params.mcmcpath << endl;
  }
  habitat.initialize(params.datapath);
  habitat.dlmwrite(params.mcmcpath);
  graph.initialize(params.datapath,habitat,params.nDemes,params.nIndiv);
  graph.dlmwrite(params.mcmcpath);
  o = graph.get_num_obsrv_demes();
  d = graph.get_num_total_demes();
  n = params.nIndiv;
  p = params.nSites;
  initialize_diffs();
  if (params.diploid) {
    Binvconst = 1.0;
    qconst = 2.0;
  } else {
    Binvconst = 4.0;
    qconst = 1.0;
  }
}
EEMS::~EEMS( ) { }
///////////////////////////////////////////
// Randraw:
void EEMS::runif_habitat(MatrixXd &Seeds) {
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
    Effcts(i) = draw.rtnorm(0.0,rateS2,-1.0*HalfInterval,HalfInterval);
  }
}
// Diffs:
void EEMS::initialize_diffs( ) {
  cerr << "[Diffs::initialize]" << endl;
  ndiv2 = (double)n/2.0; nmin1 = n-1; logn = log(n);
  J = MatrixXd::Zero(n,o);
  cvec = VectorXd::Zero(o);
  cinv = VectorXd::Zero(o);
  Diffs = MatrixXd::Zero(n,n);
  for ( int i = 0 ; i < n ; i ++ ) {
    J(i,graph.get_deme_of_indiv(i)) = 1;
    cvec(graph.get_deme_of_indiv(i)) += 1;
  }
  cinv = pow(cvec.array(),-1.0).matrix();  // cinv is the vector of inverse counts
  cmin1 = cvec; cmin1.array() -= 1;        // cmin1 is the vector of counts - 1
  int read = readMatrixXd(params.datapath + ".diffs",Diffs);
  if (read!=n*n) {
    cerr << "  Error reading dissimilarities matrix " << params.datapath + ".diffs" << endl
	 << "  Expect a " << n << "x" << n << " matrix of pairwise differences" << endl; exit(1);
  }
  cerr << "  Read dissimilarities matrix from " << params.datapath + ".diffs" << endl;
  // Check that Diffs is symmetric and full-rank
  FullPivLU<MatrixXd> lu(Diffs);
  if (lu.rank()!=n) {
    cerr << "  The dissimilarity matrix is rank-deficient" << endl; exit(1);
  }
  MatrixXd L0 = -1.0*MatrixXd::Ones(nmin1,1);
  MatrixXd L1 = MatrixXd::Identity(nmin1,nmin1);
  L = MatrixXd::Zero(nmin1,n); L << L0, L1;
  JtDhatJ = MatrixXd::Zero(o,o);
  JtDobsJ = J.transpose()*Diffs*J;
  ldLLt = logdet(L*L.transpose());
  ldLDLt = logdet(-L*Diffs*L.transpose());
  ldDiQ = ldLLt - ldLDLt;
  cerr << "[Diffs::initialize] Done." << endl << endl;
}
void EEMS::initialize(const MCMC& mcmc) {
  cerr << "[EEMS::initialize]" << endl;
  df = n;
  params.dfProposalS2 = sqrt(p);
  params.dflob = n - 1.0;
  params.dfmin = n - 1.0;
  params.dfmax = n + 1.0;
  params.dfupb = p;
  // Initialize the two Voronoi tessellations
  s2loc = draw.rinvgam(3.0,1.0);
  qtiles = draw.rnegbin(10,0.66667);
  mtiles = draw.rnegbin(10,0.66667);
  cerr << "  EEMS started with " << qtiles << " qtiles and " << mtiles << " mtiles" << endl;
  // Draw the Voronoi centers Coord uniformly within the habitat
  qSeeds = MatrixXd::Zero(qtiles,2); runif_habitat(qSeeds);
  mSeeds = MatrixXd::Zero(mtiles,2); runif_habitat(mSeeds);
  mrateS2 = draw.rinvgam(0.5,0.5);
  qrateS2 = draw.rinvgam(0.5,0.5);
  // Assign migration rates to the Voronoi tiles
  mrateMu = params.mrateMuHalfInterval*(2.0*draw.runif() - 1.0);
  // The deviation of move proposals is scaled by the habitat range
  params.mSeedsProposalS2x = params.mSeedsProposalS2 * habitat.get_xspan();
  params.mSeedsProposalS2y = params.mSeedsProposalS2 * habitat.get_yspan();
  params.qSeedsProposalS2x = params.qSeedsProposalS2 * habitat.get_xspan();
  params.qSeedsProposalS2y = params.qSeedsProposalS2 * habitat.get_yspan();
  // Assign rates to the Voronoi tiles
  qEffcts = VectorXd::Zero(qtiles); rnorm_effects(params.qEffctHalfInterval,qrateS2,qEffcts);
  mEffcts = VectorXd::Zero(mtiles); rnorm_effects(params.mEffctHalfInterval,mrateS2,mEffcts);
  // MCMC draws are stored in memory, rather than saved to disk, so it is important to thin
  niters = mcmc.num_iters_to_save();
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
  cerr << "[EEMS::initialize] Done." << endl << endl;
}
int EEMS::num_qtiles( ) const { return (qtiles); }
int EEMS::num_mtiles( ) const { return (mtiles); }
double EEMS::eval_prior( ) {
  // The parameters should always be in range
  bool inrange = true;
  for ( int i = 0 ; i < qtiles ; i++ ) {
    if (!habitat.in_point(qSeeds(i,0),qSeeds(i,1))) { inrange = false; }
  }
  for ( int i = 0 ; i < mtiles ; i++ ) {
    if (!habitat.in_point(mSeeds(i,0),mSeeds(i,1))) { inrange = false; }
  }
  if (qEffcts.cwiseAbs().minCoeff()>params.qEffctHalfInterval) { inrange = false; }
  if (mEffcts.cwiseAbs().minCoeff()>params.mEffctHalfInterval) { inrange = false; }
  if (abs(mrateMu)>params.mrateMuHalfInterval) { inrange = false; }
  if ((df<params.dfmin)||(df>params.dfmax)) { inrange = false; }
  if (!inrange) { return (-Inf); }
  double qrateS22 = 2.0*qrateS2;
  double mrateS22 = 2.0*mrateS2;
  nowpi = - log(df);
  nowpi += lgamma(params.negBiSize+mtiles) - lgamma(mtiles+1.0) + mtiles*log(params.negBiProb);
  nowpi += lgamma(params.negBiSize+qtiles) - lgamma(qtiles+1.0) + qtiles*log(params.negBiProb);
  nowpi -= (params.mrateShape/2.0+1.0)*log(mrateS22) + (params.mrateScale/2.0)/mrateS22;
  nowpi -= (params.qrateShape/2.0+1.0)*log(qrateS22) + (params.qrateScale/2.0)/qrateS22;
  nowpi -= (mtiles/2.0)*log(pi*mrateS22) + mEffcts.squaredNorm()/mrateS22;
  nowpi -= (qtiles/2.0)*log(pi*qrateS22) + qEffcts.squaredNorm()/qrateS22;
  nowpi -= (params.s2locShape/2.0+1.0)*log(s2loc) + (params.s2locScale/2.0)/s2loc;
  return (nowpi);
}
double EEMS::test_prior( ) const {
  // The parameters should always be in range
  bool inrange = true;
  for ( int i = 0 ; i < qtiles ; i++ ) {
    if (!habitat.in_point(qSeeds(i,0),qSeeds(i,1))) { inrange = false; }
  }
  for ( int i = 0 ; i < mtiles ; i++ ) {
    if (!habitat.in_point(mSeeds(i,0),mSeeds(i,1))) { inrange = false; }
  }
  if (qEffcts.cwiseAbs().minCoeff()>params.qEffctHalfInterval) { inrange = false; }
  if (mEffcts.cwiseAbs().minCoeff()>params.mEffctHalfInterval) { inrange = false; }
  if (abs(mrateMu)>params.mrateMuHalfInterval) { inrange = false; }
  if ((df<params.dfmin)||(df>params.dfmax)) { inrange = false; }
  if (!inrange) { return (-Inf); }
  double qrateS22 = 2.0*qrateS2;
  double mrateS22 = 2.0*mrateS2;
  double logpi = - log(df);
  logpi += lgamma(params.negBiSize+mtiles) - lgamma(mtiles+1.0) + mtiles*log(params.negBiProb);
  logpi += lgamma(params.negBiSize+qtiles) - lgamma(qtiles+1.0) + qtiles*log(params.negBiProb);
  logpi -= (params.mrateShape/2.0+1.0)*log(mrateS22) + (params.mrateScale/2.0)/mrateS22;
  logpi -= (params.qrateShape/2.0+1.0)*log(qrateS22) + (params.qrateScale/2.0)/qrateS22;
  logpi -= (mtiles/2.0)*log(pi*mrateS22) + mEffcts.squaredNorm()/mrateS22;
  logpi -= (qtiles/2.0)*log(pi*qrateS22) + qEffcts.squaredNorm()/qrateS22;
  logpi -= (params.s2locShape/2.0+1.0)*log(s2loc) + (params.s2locScale/2.0)/s2loc;
  return (logpi);
}
double EEMS::eval_likelihood( ) {
  graph.pdist2(mSeeds,nowmColors);
  graph.pdist2(qSeeds,nowqColors);
  calc_q(nowqColors,qEffcts,nowq);
  calc_Binv(nowmColors,mEffcts,mrateMu,nowBinv);
  double trDinvQxD, ll_partdf;
  nowll = EEMS_wishpdfln(nowBinv,nowq,s2loc,df,trDinvQxD,ll_partdf);
  nowtrDinvQxD = trDinvQxD;
  nowll_partdf = ll_partdf;
  return (nowll);
}
double EEMS::test_likelihood( ) const {
  VectorXi mColors, qColors;
  graph.pdist2(mSeeds,mColors);
  graph.pdist2(qSeeds,qColors);
  VectorXd q = VectorXd::Zero(d);
  for ( int alpha = 0 ; alpha < d ; alpha++ ) {
    double log10m = qEffcts(qColors(alpha)); // qrateMu = 0.0
    q(alpha) = pow(10.0,log10m);
  }
  MatrixXd M = MatrixXd::Zero(d,d);
  int alpha, beta;
  for ( int edge = 0 ; edge < graph.get_num_edges() ; edge++ ) {
    graph.get_edge(edge,alpha,beta);
    double log10m1 = mEffcts(mColors(alpha)) + mrateMu;
    double log10m2 = mEffcts(mColors(beta)) + mrateMu;
    M(alpha,beta) = 0.5 * pow(10.0,log10m1) + 0.5 * pow(10.0,log10m2);
  }
  MatrixXd Delta = expected_dissimilarities(J, (M + M.transpose()) * Binvconst, q * qconst);
  double logll = wishpdfln(-L*Diffs*L.transpose(),
			   -L*Delta*L.transpose()*s2loc/df,df);
  return (logll);
}
void EEMS::calc_q(const VectorXi &qColors0, const VectorXd &qEffcts0, VectorXd &q0) const {
  if (q0.size()!=o) { q0.resize(o); }
  for ( int alpha = 0 ; alpha < o ; alpha++ ) {
    double log10m = qEffcts0(qColors0(alpha));
    q0(alpha) = pow(10.0,log10m);
  }
  q0 *= qconst;
}
void EEMS::calc_Binv(const VectorXi &mColors0, const VectorXd &mEffcts0, const double mrateMu0, MatrixXd &Binv0) const {
  if (Binv0.rows()!=n||Binv0.cols()!=d) { Binv0.resize(d,d); }
  vector<Tri> coefficients;
  int alpha, beta;
  for ( int edge = 0 ; edge < graph.get_num_edges() ; edge++ ) {
    graph.get_edge(edge,alpha,beta);
    double log10m1 = mrateMu0 + mEffcts0(mColors0(alpha));
    double log10m2 = mrateMu0 + mEffcts0(mColors0(beta));
    double m12 = 0.5 * pow(10.0,log10m1) + 0.5 * pow(10.0,log10m2);
    coefficients.push_back(Tri(alpha,beta,m12));
    coefficients.push_back(Tri(beta,alpha,m12));
  }
  SpMat sparseM(d,d);
  sparseM.setFromTriplets(coefficients.begin(),coefficients.end());
  MatrixXd M = MatrixXd(sparseM);
  MatrixXd Hinv = - M; Hinv.diagonal() += M.rowwise().sum(); Hinv.array() += 1.0;
  if (o==d) {
    Binv0 = -0.5 * Hinv;
  } else {
    Binv0 = -0.5 * Hinv.topLeftCorner(o,o);
    Binv0 += 0.5 * Hinv.topRightCorner(o,d-o) *
      Hinv.bottomRightCorner(d-o,d-o).selfadjointView<Lower>().llt().solve(Hinv.bottomLeftCorner(d-o,o));
  }
  Binv0 *= Binvconst;
}
double EEMS::EEMS_wishpdfln(const MatrixXd &Binv, const VectorXd &q, const double s2loc, const double df,
			    double &trDinvQxD, double &ll_partdf) const {
  double df2 = 0.5 * df;
  MatrixXd T = Binv.selfadjointView<Lower>().ldlt().solve(MatrixXd::Identity(o,o));
  T *= cvec.asDiagonal(); T -= q.asDiagonal(); // Now T = B*C - W
  PartialPivLU<MatrixXd> lu(T);
  VectorXd qinv = pow(q.array(),-1.0);
  T.noalias() = T*qinv.asDiagonal()*cinv.asDiagonal();
  T += cinv.asDiagonal();                   // Now T = B*Winv
  MatrixXd X = lu.solve(T);
  VectorXd Xc_qinv = X*cvec - qinv;
  double oDinvo = cvec.dot(Xc_qinv);
  double oDiDDi = Xc_qinv.transpose()*JtDobsJ*Xc_qinv;
  trDinvQxD = trace_AxB(X,JtDobsJ) - oDiDDi/oDinvo;
  ll_partdf = logn - log(abs(oDinvo))
    + cmin1.dot(qinv.array().log().matrix())
    - lu.matrixLU().diagonal().array().abs().log().sum()
    - trDinvQxD/s2loc - nmin1*log(s2loc) - ldDiQ;
  return (df2 * ll_partdf + nmin1*df2*log(df2) - mvgammaln(df2,nmin1) - ndiv2*ldLDLt);
}
double EEMS::eval_proposal_qEffcts(Proposal &proposal) const {
  VectorXd newqEffcts = qEffcts;
  newqEffcts(proposal.qTile) = proposal.newqEffct;
  calc_q(nowqColors,newqEffcts,proposal.newq);
  return (EEMS_wishpdfln(nowBinv,proposal.newq,s2loc,df,proposal.newtrDinvQxD,proposal.newll_partdf));
}
double EEMS::eval_proposal_qSeeds(Proposal &proposal) const {
  MatrixXd newqSeeds = qSeeds;
  newqSeeds(proposal.qTile,0) = proposal.newqSeedx;
  newqSeeds(proposal.qTile,1) = proposal.newqSeedy;
  graph.pdist2(newqSeeds,proposal.newqColors);
  calc_q(proposal.newqColors,qEffcts,proposal.newq);
  return (EEMS_wishpdfln(nowBinv,proposal.newq,s2loc,df,proposal.newtrDinvQxD,proposal.newll_partdf));
}
double EEMS::eval_birthdeath_qVoronoi(Proposal &proposal) const {
  graph.pdist2(proposal.newqSeeds,proposal.newqColors);
  calc_q(proposal.newqColors,proposal.newqEffcts,proposal.newq);
  return (EEMS_wishpdfln(nowBinv,proposal.newq,s2loc,df,proposal.newtrDinvQxD,proposal.newll_partdf));
}
double EEMS::eval_proposal_mEffcts(Proposal &proposal) const {
  VectorXd newmEffcts = mEffcts;
  newmEffcts(proposal.mTile) = proposal.newmEffct;
  calc_Binv(nowmColors,newmEffcts,mrateMu,proposal.newBinv);
  return (EEMS_wishpdfln(proposal.newBinv,nowq,s2loc,df,proposal.newtrDinvQxD,proposal.newll_partdf));
}
double EEMS::eval_proposal_mrateMu(Proposal &proposal) const {
  calc_Binv(nowmColors,mEffcts,proposal.newmrateMu,proposal.newBinv);
  return (EEMS_wishpdfln(proposal.newBinv,nowq,s2loc,df,proposal.newtrDinvQxD,proposal.newll_partdf));
}
double EEMS::eval_proposal_mSeeds(Proposal &proposal) const {
  MatrixXd newmSeeds = mSeeds;
  newmSeeds(proposal.mTile,0) = proposal.newmSeedx;
  newmSeeds(proposal.mTile,1) = proposal.newmSeedy;
  graph.pdist2(newmSeeds,proposal.newmColors);
  calc_Binv(proposal.newmColors,mEffcts,mrateMu,proposal.newBinv);
  return (EEMS_wishpdfln(proposal.newBinv,nowq,s2loc,df,proposal.newtrDinvQxD,proposal.newll_partdf));
}
double EEMS::eval_birthdeath_mVoronoi(Proposal &proposal) const {
  graph.pdist2(proposal.newmSeeds,proposal.newmColors);
  calc_Binv(proposal.newmColors,proposal.newmEffcts,mrateMu,proposal.newBinv);
  return (EEMS_wishpdfln(proposal.newBinv,nowq,s2loc,df,proposal.newtrDinvQxD,proposal.newll_partdf));
}
///////////////////////////////////////////
// Proposals:
void EEMS::update_s2loc( ) {
  double df2 = 0.5 * df;
  nowll_partdf += nowtrDinvQxD/s2loc + nmin1*log(s2loc);
  nowpi += (params.s2locShape/2.0+1.0)*log(s2loc) + (params.s2locScale/2.0)/s2loc;
  s2loc = draw.rinvgam( 0.5 * ( params.s2locShape + df*nmin1 ),
			0.5 * ( params.s2locScale + df*nowtrDinvQxD ) );
  nowll_partdf -= nowtrDinvQxD/s2loc + nmin1*log(s2loc);
  nowpi -= (params.s2locShape/2.0+1.0)*log(s2loc) + (params.s2locScale/2.0)/s2loc;
  nowll = df2 * nowll_partdf + nmin1*df2*log(df2) - mvgammaln(df2,nmin1) - ndiv2*ldLDLt;
}
void EEMS::propose_df(Proposal &proposal) {
  double newdf = draw.rnorm(df,params.dfProposalS2);
  double newdf2 = 0.5 * newdf;
  proposal.type = 7;
  proposal.newdf = newdf;
  if ( (newdf>params.dfmin) && (newdf<params.dfmax) ) {
    proposal.newpi = nowpi + log(df) - log(newdf);
    proposal.newll = newdf2 * nowll_partdf + nmin1*newdf2*log(newdf2) - mvgammaln(newdf2,nmin1) - ndiv2*ldLDLt;
  } else {
    proposal.newpi = -Inf;
    proposal.newll = -Inf;
  }
}
void EEMS::propose_qEffcts(Proposal &proposal, const MCMC &mcmc) {
  double curqEffct = qEffcts(mcmc.qTile);
  double newqEffct = draw.rnorm(curqEffct,params.qEffctProposalS2);
  proposal.type = 0;
  proposal.qTile = mcmc.qTile;
  proposal.newqEffct = newqEffct;
  if ( abs(newqEffct) < params.qEffctHalfInterval ) {
    proposal.newpi = nowpi - (newqEffct * newqEffct - curqEffct * curqEffct) / (2.0 * qrateS2);
    proposal.newll = eval_proposal_qEffcts(proposal);
  } else {
    proposal.newpi = -Inf;
    proposal.newll = -Inf;
  }
}
void EEMS::propose_mEffcts(Proposal &proposal, const MCMC &mcmc) {
  double curmEffct = mEffcts(mcmc.mTile);
  double newmEffct = draw.rnorm(curmEffct,params.mEffctProposalS2);
  proposal.type = 3;
  proposal.mTile = mcmc.mTile;
  proposal.newmEffct = newmEffct;
  if ( abs(newmEffct) < params.mEffctHalfInterval ) {
    proposal.newpi = nowpi - (newmEffct * newmEffct - curmEffct * curmEffct) / (2.0 * mrateS2);
    proposal.newll = eval_proposal_mEffcts(proposal);
  } else {
    proposal.newpi = -Inf;
    proposal.newll = -Inf;
  }
}
void EEMS::propose_mrateMu(Proposal &proposal) {
  double newmrateMu = draw.rnorm(mrateMu,params.mrateMuProposalS2);
  proposal.type = 4;
  proposal.newmrateMu = newmrateMu;
  if ( abs(newmrateMu) < params.mrateMuHalfInterval ) {
    proposal.newpi = nowpi;
    proposal.newll = eval_proposal_mrateMu(proposal);
  } else {
    proposal.newpi = -Inf;
    proposal.newll = -Inf;
  }
}
void EEMS::move_qVoronoi(Proposal &proposal, const MCMC &mcmc) {
  double newqSeedx = draw.rnorm(qSeeds(mcmc.qTile,0),params.qSeedsProposalS2x);
  double newqSeedy = draw.rnorm(qSeeds(mcmc.qTile,1),params.qSeedsProposalS2y);
  proposal.type = 1;
  proposal.qTile = mcmc.qTile;
  proposal.newqSeedx = newqSeedx;
  proposal.newqSeedy = newqSeedy;
  if (habitat.in_point(newqSeedx,newqSeedy)) {
    proposal.newpi = nowpi;
    proposal.newll = eval_proposal_qSeeds(proposal);
  } else {
    proposal.newpi = -Inf;
    proposal.newll = -Inf;
  }
}
void EEMS::move_mVoronoi(Proposal &proposal, const MCMC &mcmc) {
  double newmSeedx = draw.rnorm(mSeeds(mcmc.mTile,0),params.mSeedsProposalS2x);
  double newmSeedy = draw.rnorm(mSeeds(mcmc.mTile,1),params.mSeedsProposalS2y);
  proposal.type = 5;
  proposal.mTile = mcmc.mTile;
  proposal.newmSeedx = newmSeedx;
  proposal.newmSeedy = newmSeedy;
  if (habitat.in_point(newmSeedx,newmSeedy)) {
    proposal.newpi = nowpi;
    proposal.newll = eval_proposal_mSeeds(proposal);
  } else {
    proposal.newpi = -Inf;
    proposal.newll = -Inf;
  }
}
void EEMS::birthdeath_qVoronoi(Proposal &proposal, const MCMC &mcmc) {
  int newqtiles = qtiles;
  double u = draw.runif();
  double pBirth = 0.5;
  double pDeath = 0.5;
  proposal.newqEffcts = qEffcts;
  proposal.newqSeeds = qSeeds;
  // Since there must be at least one tile, rule out a death proposal
  if (qtiles==1) { u = 0.0; }
  if (u < 0.5) {
    // Propose birth
    newqtiles++;
    if (qtiles==1) { pBirth = 1.0; }
    MatrixXd newqSeed = MatrixXd::Zero(1,2);
    VectorXd newqEffct = VectorXd::Zero(1);
    runif_habitat(newqSeed);
    rnorm_effects(params.qEffctHalfInterval,qrateS2,newqEffct);
    insertRow(proposal.newqSeeds,newqSeed.row(0));
    insertElem(proposal.newqEffcts,newqEffct(0));
    proposal.newpi = nowpi + log(pDeath/pBirth) + log((qtiles+params.negBiSize)/(newqtiles/params.negBiProb));
  } else {
    // Propose death
    newqtiles--;
    if (qtiles==2) { pBirth = 1.0; }
    int qtileToRemove = draw.runif_int(0,newqtiles);
    removeRow(proposal.newqSeeds,qtileToRemove);
    removeElem(proposal.newqEffcts,qtileToRemove);
    proposal.newpi = nowpi + log(pBirth/pDeath) + log((qtiles/params.negBiProb)/(newqtiles+params.negBiSize));
  }
  proposal.type = 2;
  proposal.newqtiles = newqtiles;
  proposal.newll = eval_birthdeath_qVoronoi(proposal);
}
void EEMS::birthdeath_mVoronoi(Proposal &proposal, const MCMC &mcmc) {
  int newmtiles = mtiles;
  double u = draw.runif();
  double pBirth = 0.5;
  double pDeath = 0.5;
  proposal.newmEffcts = mEffcts;
  proposal.newmSeeds = mSeeds;
  // Since there must be at least one tile, rule out a death proposal
  if (mtiles==1) { u = 0.0; }
  if (u < 0.5) {
    // Propose birth
    newmtiles++;
    if (mtiles==1) { pBirth = 1.0; }
    MatrixXd newmSeed = MatrixXd::Zero(1,2);
    VectorXd newmEffct = VectorXd::Zero(1);
    runif_habitat(newmSeed);
    rnorm_effects(params.mEffctHalfInterval,mrateS2,newmEffct);
    insertRow(proposal.newmSeeds,newmSeed.row(0));
    insertElem(proposal.newmEffcts,newmEffct(0));
    proposal.newpi = nowpi + log(pDeath/pBirth) + log((mtiles+params.negBiSize)/(newmtiles/params.negBiProb));
  } else {
    // Propose death
    newmtiles--;
    if (mtiles==2) { pBirth = 1.0; }
    int mtileToRemove = draw.runif_int(0,newmtiles);
    removeRow(proposal.newmSeeds,mtileToRemove);
    removeElem(proposal.newmEffcts,mtileToRemove);
    proposal.newpi = nowpi + log(pBirth/pDeath) + log((mtiles/params.negBiProb)/(newmtiles+params.negBiSize));
  }
  proposal.type = 6;
  proposal.newmtiles = newmtiles;
  proposal.newll = eval_birthdeath_mVoronoi(proposal);
}
void EEMS::update_df_support(const MCMC &mcmc) {
  if (mcmc.currIter > (mcmc.numBurnIter/3)) {
    params.dfmax = params.dfupb;
    params.dfmin = params.dflob;
  }
}
void EEMS::update_hyperparams( ) {
  double qrateS22 = 2.0*qrateS2;
  double mrateS22 = 2.0*mrateS2;
  double SSq = qEffcts.squaredNorm();
  double SSm = mEffcts.squaredNorm();
  nowpi += (params.mrateShape/2.0+1.0)*log(mrateS22) + (params.mrateScale/2.0)/mrateS22;
  nowpi += (params.qrateShape/2.0+1.0)*log(qrateS22) + (params.qrateScale/2.0)/qrateS22;
  nowpi += (mtiles/2.0)*log(pi*mrateS22) + SSm/mrateS22;
  nowpi += (qtiles/2.0)*log(pi*qrateS22) + SSq/qrateS22;
  qrateS2 = draw.rinvgam(0.5 * (params.qrateShape + qtiles), 0.5 * (params.qrateScale + SSq));
  mrateS2 = draw.rinvgam(0.5 * (params.mrateShape + mtiles), 0.5 * (params.mrateScale + SSm));
  qrateS22 = 2.0*qrateS2;
  mrateS22 = 2.0*mrateS2;
  nowpi -= (params.mrateShape/2.0+1.0)*log(mrateS22) + (params.mrateScale/2.0)/mrateS22;
  nowpi -= (params.qrateShape/2.0+1.0)*log(qrateS22) + (params.qrateScale/2.0)/qrateS22;
  nowpi -= (mtiles/2.0)*log(pi*mrateS22) + SSm/mrateS22;
  nowpi -= (qtiles/2.0)*log(pi*qrateS22) + SSq/qrateS22;
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
  double alpha = min(0.0,(proposal.newpi-nowpi) + (proposal.newll-nowll));
  int qTile, mTile;
  if ( log(u) < alpha ) {
    // Update parameters:
    switch (proposal.type) {
    case 0:
      qTile = proposal.qTile;
      qEffcts(qTile) = proposal.newqEffct;
      nowq = proposal.newq;
      break;
    case 1:
      qTile = proposal.qTile;
      qSeeds(qTile,0) = proposal.newqSeedx;
      qSeeds(qTile,1) = proposal.newqSeedy;
      nowq = proposal.newq;
      nowqColors = proposal.newqColors;
      break;
    case 2:
      qSeeds = proposal.newqSeeds;
      qtiles = proposal.newqtiles;
      nowq = proposal.newq;
      qEffcts = proposal.newqEffcts;
      nowqColors = proposal.newqColors;
      break;
    case 3:
      mTile = proposal.mTile;
      mEffcts(mTile) = proposal.newmEffct;
      nowBinv = proposal.newBinv;
      break;
    case 4:
      mrateMu = proposal.newmrateMu;
      nowBinv = proposal.newBinv;
      break;
    case 5:
      mTile = proposal.mTile;
      mSeeds(mTile,0) = proposal.newmSeedx;
      mSeeds(mTile,1) = proposal.newmSeedy;
      nowBinv = proposal.newBinv;
      nowmColors = proposal.newmColors;
      break;
    case 6:
      mSeeds = proposal.newmSeeds;
      mtiles = proposal.newmtiles;
      mEffcts = proposal.newmEffcts;
      nowBinv = proposal.newBinv;
      nowmColors = proposal.newmColors;
      break;
    default:
      df = proposal.newdf;
    }
    nowpi = proposal.newpi;
    nowll = proposal.newll;
    if (proposal.type!=7) { // 7 means a df proposal
      nowtrDinvQxD = proposal.newtrDinvQxD;
      nowll_partdf = proposal.newll_partdf;
    }
    if (proposal.type==2||proposal.type==6) {
      nowpi = eval_prior();
    }
    return true;
  }
  proposal.newpi = nowpi;
  proposal.newll = nowll;
  return false;
}
///////////////////////////////////////////
// Save results:
bool EEMS::save_iteration(const int iter) {
  if ((iter<0)||(iter>=niters)) { return false; }
  mcmcthetas(iter,0) = s2loc;
  mcmcthetas(iter,1) = df;
  mcmcqhyper(iter,0) = 0.0;
  mcmcqhyper(iter,1) = qrateS2;
  mcmcmhyper(iter,0) = mrateMu;
  mcmcmhyper(iter,1) = mrateS2;
  mcmcpilogl(iter,0) = nowpi;
  mcmcpilogl(iter,1) = nowll;
  mcmcqtiles(iter) = qtiles;
  mcmcmtiles(iter) = mtiles;
  for ( int t = 0 ; t < qtiles ; t++ ) {
    mcmcqRates.push_back(pow(10.0,qEffcts(t)));
  }
  for ( int t = 0 ; t < qtiles ; t++ ) {
    mcmcwCoord.push_back(qSeeds(t,0));
  }
  for ( int t = 0 ; t < qtiles ; t++ ) {
    mcmczCoord.push_back(qSeeds(t,1));
  }
  for ( int t = 0 ; t < mtiles ; t++ ) {
    mcmcmRates.push_back(pow(10.0,mEffcts(t) + mrateMu));
  }
  for ( int t = 0 ; t < mtiles ; t++ ) {
    mcmcxCoord.push_back(mSeeds(t,0));
  }
  for ( int t = 0 ; t < mtiles ; t++ ) {
    mcmcyCoord.push_back(mSeeds(t,1));
  }
  MatrixXd B = nowBinv.inverse();
  VectorXd h = B.diagonal();
  B -= 0.5 * h.replicate(1,o);
  B -= 0.5 * h.transpose().replicate(o,1);
  B += 0.5 * nowq.replicate(1,o);
  B += 0.5 * nowq.transpose().replicate(o,1);
  JtDhatJ += s2loc * B;
  return true;
}
bool EEMS::output_results(const MCMC &mcmc) const {
  ofstream out; bool err;
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
  err = dlmcell(params.mcmcpath + "/mcmcmrates.txt",mcmcmtiles,mcmcmRates); if (err) { return false; }
  err = dlmcell(params.mcmcpath + "/mcmcxcoord.txt",mcmcmtiles,mcmcxCoord); if (err) { return false; }
  err = dlmcell(params.mcmcpath + "/mcmcycoord.txt",mcmcmtiles,mcmcyCoord); if (err) { return false; }
  err = dlmcell(params.mcmcpath + "/mcmcqrates.txt",mcmcqtiles,mcmcqRates); if (err) { return false; }
  err = dlmcell(params.mcmcpath + "/mcmcwcoord.txt",mcmcqtiles,mcmcwCoord); if (err) { return false; }
  err = dlmcell(params.mcmcpath + "/mcmczcoord.txt",mcmcqtiles,mcmczCoord); if (err) { return false; }
  out.open((params.mcmcpath + "/eemsrun.txt").c_str(),ofstream::out);
  if (!out.is_open( )) { return false; }
  out << "[EEMS::Params] Random seed = "
      << draw.get_seed( ) << endl << endl;
  params.dlmwrite(out);
  out << "Acceptance proportions:" << endl;
  mcmc.output_proportions(out); out << endl;
  out << "Final log prior = " << nowpi << endl
      << "Final log llike = " << nowll << endl;
  out.close( );
  return true;
}
void EEMS::report_iteration(const int iter) const {
  cerr << fixed << setprecision(2)
       << "and effective degrees of freedom = " << df << endl
       << "        number of qVoronoi tiles = " << qtiles << endl
       << "        number of mVoronoi tiles = " << mtiles << endl
       << "          Log prior = " << nowpi << endl
       << "          Log llike = " << nowll << endl;
}
void EEMS::check_ll_computation( ) const {
  double pi0 = test_prior( );
  double ll0 = test_likelihood( );
  if ((abs(nowpi-pi0)/abs(pi0)>1e-12)||
      (abs(nowll-ll0)/abs(ll0)>1e-12)) {
    cerr << "[EEMS::testing]   |ll0-ll|/|ll0| = " << abs(nowll - ll0)/abs(ll0) << endl;
    cerr << "[EEMS::testing]   |pi0-pi|/|pi0| = " << abs(nowpi - pi0)/abs(pi0) << endl;
    exit(1);
  }
}
