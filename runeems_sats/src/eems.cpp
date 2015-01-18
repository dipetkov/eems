
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
    Bconst = 1.0;
    qconst = 2.0;
  } else {
    Bconst = 0.25;
    qconst = 1.0;
  }
}
EEMS::~EEMS( ) { }
///////////////////////////////////////////
// Randraw:
double EEMS::runif( ) { return (draw.runif( )); }
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
  int alleles_to_read = 0;
  if (params.diploid) {
    alleles_to_read = 2*p;
  } else {
    alleles_to_read = p;
  }
  MatrixXd Sites = readMatrixXd(params.datapath + ".sites");
  if ((Sites.rows()!=n)||(Sites.cols()!=alleles_to_read)) {
    cerr << "  Error reading genotype data matrix " << params.datapath + ".sites" << endl
	 << "  Expect a " << n << "x" << alleles_to_read << " matrix of allele copies" << endl;
    exit(1);
  }
  cerr << "  Read genotype data from " << params.datapath + ".diffs" << endl;
  ///////////////////////////////////////
  J = MatrixXd::Zero(n,o);
  for ( int i = 0 ; i < n ; i ++ ) {
    J(i,graph.get_deme_of_indiv(i)) = 1;
  }
  MatrixXd Diffs_allSites = MatrixXd::Zero(n,n);
  MatrixXd Pairs_allSites = MatrixXd::Zero(n,n);
  ///////////////////////////////////////
  logn.resize(p);
  nmin1.resize(p);
  ll_partdf = 0.0;
  for ( int i = 0 ; i < p ; i ++ ) {
    VectorXd z = VectorXd::Zero(n);
    if (params.diploid) {
      VectorXd a1 = Sites.col(2*i);
      VectorXd a2 = Sites.col(2*i+1);
      if (((a1.array()<0) != (a2.array()<0)).maxCoeff()) {
	cerr << "  Error processing genotypes: one allele is missing but the other is not." << endl;
	exit(1);
      }
      z = 0.5 * a1 + 0.5 * a2;
    } else {
      z = Sites.col(i);
    }
    int ni = (z.array()>=0).count();
    // Copy into a vector with the unobserved samples removed.
    int j = 0;
    VectorXd zi = VectorXd::Zero(ni);
    VectorXi oi = VectorXi::Zero(ni);
    for( int j2 = 0 ; j2 < ni ; j2++ ) {
      while (z(j)<0) { j++; }
      zi(j2) = z(j);
      oi(j2) = j++;
    }
    // Indicates whether the sample's allele is observed or missing
    VectorXd I = (z.array()<0).select(VectorXd::Zero(n),VectorXd::Ones(n));
    MatrixXd Pairs = I * I.transpose();
    MatrixXd Ji = slice(J,oi,VectorXi::LinSpaced(o,0,o-1));
    MatrixXd Li = MatrixXd::Zero(ni-1,ni);
    Li << - MatrixXd::Ones(ni-1,1), MatrixXd::Identity(ni-1,ni-1);
    Z.push_back(zi);
    O.push_back(oi);
    L.push_back(Li);
    // Similarity matrix (with rank 1)
    MatrixXd Si = zi*zi.transpose();
    // Dissimilarity matrix
    MatrixXd Di = - 2.0 * Si;
    Di += Si.diagonal().replicate(1,ni);
    Di += Si.diagonal().transpose().replicate(ni,1);
    ll_partdf += logdet(Li*Li.transpose()) + (ni-1.0) * ( pseudologdet( - Li*Di*Li.transpose(), 1) + log_2 + log_pi );
    logn(i) = log(ni);
    nmin1(i) = ni - 1;
    VectorXd ci = Ji.colwise().sum();
    int o0 = (ci.array()>0).count();
    VectorXi si = VectorXi::Zero(o0);
    j = 0;
    for (int j2 = 0 ; j2 < o0 ; j2 ++ ) {
      while (ci(j)==0) { j++; }
      si(j2) = j++;
    }
    ci = slice(ci,si);
    cvec.push_back(ci);
    cinv.push_back(pow(ci.array(), - 1.0));
    cmin1.push_back(ci.array() - 1.0);
    Diffs.push_back(Di);
    MatrixXd JtDJ = Ji.transpose()*Di*Ji;
    JtDobsJ.push_back(slice(JtDJ,si,si));
    ovec.push_back(si);
    /////////////////////////////////////
    double mu = zi.mean();
    zi.array() -=  mu;
    double sd = sqrt( zi.squaredNorm() / (ni - 1.0) );
    z.array() -= mu;
    z.array() /= sd;
    Si = z*z.transpose();
    // Dissimilarity matrix
    Di = - 2.0 * Si;
    Di += Si.diagonal().replicate(1,n);
    Di += Si.diagonal().transpose().replicate(n,1);
    // Set rows/cols that correspond to missing samples to 0s
    Di.array() *= Pairs.array();
    Diffs_allSites += Di;
    Pairs_allSites += Pairs;
  }
  Diffs_allSites.array() /= Pairs_allSites.array();
  JtDobsJ_allSites = J.transpose()*Diffs_allSites*J;
  JtDhatJ_allSites = MatrixXd::Zero(o,o);
  c_allSites = J.colwise().sum();
  cerr << "[Diffs::initialize] Done." << endl << endl;
}
void EEMS::initialize(const MCMC& mcmc) {
  cerr << "[EEMS::initialize]" << endl;
  s2loc = VectorXd::Zero(p);
  for (int i = 0 ; i < p ; i++ ) {
    s2loc(i) = draw.rinvgam(3.0,1.0);
  }
  // Initialize the two Voronoi tessellations
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
  mcmcthetas = MatrixXd::Zero(niters,p);
  mcmcmhyper = MatrixXd::Zero(niters,2);
  mcmcqhyper = MatrixXd::Zero(niters,2);
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
  // The parameters should be in range by construction
  double qrateS22 = 2.0*qrateS2;
  double mrateS22 = 2.0*mrateS2;
  nowpi  = lgamma(params.negBiSize+mtiles) - lgamma(mtiles+1.0) + mtiles*log(params.negBiProb);
  nowpi += lgamma(params.negBiSize+qtiles) - lgamma(qtiles+1.0) + qtiles*log(params.negBiProb);
  nowpi -= (params.mrateShape/2.0+1.0)*log(mrateS22) + (params.mrateScale/2.0)/mrateS22;
  nowpi -= (params.qrateShape/2.0+1.0)*log(qrateS22) + (params.qrateScale/2.0)/qrateS22;
  nowpi -= (mtiles/2.0)*log(pi*mrateS22) + mEffcts.squaredNorm()/mrateS22;
  nowpi -= (qtiles/2.0)*log(pi*qrateS22) + qEffcts.squaredNorm()/qrateS22;
  nowpi -= (params.s2locShape/2.0+1.0)*s2loc.array().log().sum();
  nowpi -= (params.s2locScale/2.0)*pow(s2loc.array(),-1.0).sum();
  return (nowpi);
}
double EEMS::test_prior( ) const {
  // The parameters should be in range by construction
  double qrateS22 = 2.0*qrateS2;
  double mrateS22 = 2.0*mrateS2;
  double logpi;
  logpi  = lgamma(params.negBiSize+mtiles) - lgamma(mtiles+1.0) + mtiles*log(params.negBiProb);
  logpi += lgamma(params.negBiSize+qtiles) - lgamma(qtiles+1.0) + qtiles*log(params.negBiProb);
  logpi -= (params.mrateShape/2.0+1.0)*log(mrateS22) + (params.mrateScale/2.0)/mrateS22;
  logpi -= (params.qrateShape/2.0+1.0)*log(qrateS22) + (params.qrateScale/2.0)/qrateS22;
  logpi -= (mtiles/2.0)*log(pi*mrateS22) + mEffcts.squaredNorm()/mrateS22;
  logpi -= (qtiles/2.0)*log(pi*qrateS22) + qEffcts.squaredNorm()/qrateS22;
  for ( int i = 0 ; i < p ; i++ ) {
    logpi -= (params.s2locShape/2.0+1.0)*log(s2loc(i)) + (params.s2locScale/2.0)/s2loc(i);
  }
  return (logpi);
}
double EEMS::eval_likelihood( ) {
  graph.pdist2(mSeeds,nowmColors);
  graph.pdist2(qSeeds,nowqColors);
  calc_q(nowqColors,qEffcts,nowq);
  calc_B(nowmColors,mEffcts,mrateMu,nowB);
  VectorXd trDinvQxD = VectorXd::Zero(p);
  nowll = EEMS_wishpdfln(nowB,nowq,s2loc,trDinvQxD);
  nowtrDinvQxD = trDinvQxD;
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
  MatrixXd Delta = expected_dissimilarities(J, (M + M.transpose()) / Bconst, q * qconst);
  double logll = 0.0;
  for ( int i = 0 ; i < p ; i++ ) {
    logll += pseudowishpdfln(-L[i]*Diffs[i]*L[i].transpose(),
			     -L[i]*slice(Delta,O[i],O[i])*L[i].transpose()*s2loc(i),1);
  }
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
void EEMS::calc_B(const VectorXi &mColors0, const VectorXd &mEffcts0, const double mrateMu0, MatrixXd &B0) const {
  if (B0.rows()!=n||B0.cols()!=d) { B0.resize(d,d); }
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
  MatrixXd Binv0;
  if (o==d) {
    Binv0 = -0.5 * Hinv;
  } else {
    Binv0 = -0.5 * Hinv.topLeftCorner(o,o);
    Binv0 += 0.5 * Hinv.topRightCorner(o,d-o) *
      Hinv.bottomRightCorner(d-o,d-o).selfadjointView<Lower>().llt().solve(Hinv.bottomLeftCorner(d-o,o));
  }
  B0 = Bconst * Binv0.inverse();
  VectorXd h = B0.diagonal();
  B0 -= 0.5 * h.replicate(1,o);
  B0 -= 0.5 * h.transpose().replicate(o,1);
}
double EEMS::EEMS_wishpdfln(const MatrixXd &B, const VectorXd &q, const VectorXd &s2loc, VectorXd &trDinvQxD) const {
  VectorXd ldetDinvQ = VectorXd::Zero(p);
  if (trDinvQxD.size() != p) { trDinvQxD.resize(p); }
  for ( int i = 0 ; i < p ; i++ ) {
    VectorXd qi = slice(q,ovec[i]);
    VectorXd qiinv = pow(qi.array(),-1.0);
    MatrixXd T = slice(B,ovec[i],ovec[i]);
    T *= cvec[i].asDiagonal(); T -= qi.asDiagonal(); // Now T = B*C - W
    PartialPivLU<MatrixXd> lu(T);
    T.noalias() = T*qiinv.asDiagonal()*cinv[i].asDiagonal();
    T += cinv[i].asDiagonal();                       // Now T = B*Winv
    MatrixXd X = lu.solve(T);
    VectorXd Xc_qinv = X*cvec[i] - qiinv;
    double oDinvo = cvec[i].dot(Xc_qinv);
    double oDiDDi = Xc_qinv.transpose()*JtDobsJ[i]*Xc_qinv;
    ldetDinvQ(i) = logn(i) - log(abs(oDinvo))
      + cmin1[i].dot(qiinv.array().log().matrix())
      - lu.matrixLU().diagonal().array().abs().log().sum();
    trDinvQxD(i) = trace_AxB(X,JtDobsJ[i]) - oDiDDi/oDinvo;
  }
  return ( 0.5 * ( ldetDinvQ.sum() - (trDinvQxD.array() / s2loc.array()).sum() -
		   (nmin1.array() * s2loc.array().log()).sum() - ll_partdf ) );
}
double EEMS::eval_proposal_qEffcts(Proposal &proposal) const {
  VectorXd newqEffcts = qEffcts;
  newqEffcts(proposal.qTile) = proposal.newqEffct;
  calc_q(nowqColors,newqEffcts,proposal.newq);
  return (EEMS_wishpdfln(nowB,proposal.newq,s2loc,proposal.newtrDinvQxD));
}
double EEMS::eval_proposal_qSeeds(Proposal &proposal) const {
  MatrixXd newqSeeds = qSeeds;
  newqSeeds(proposal.qTile,0) = proposal.newqSeedx;
  newqSeeds(proposal.qTile,1) = proposal.newqSeedy;
  graph.pdist2(newqSeeds,proposal.newqColors);
  calc_q(proposal.newqColors,qEffcts,proposal.newq);
  return (EEMS_wishpdfln(nowB,proposal.newq,s2loc,proposal.newtrDinvQxD));
}
double EEMS::eval_birthdeath_qVoronoi(Proposal &proposal) const {
  graph.pdist2(proposal.newqSeeds,proposal.newqColors);
  calc_q(proposal.newqColors,proposal.newqEffcts,proposal.newq);
  return (EEMS_wishpdfln(nowB,proposal.newq,s2loc,proposal.newtrDinvQxD));
}
double EEMS::eval_proposal_mEffcts(Proposal &proposal) const {
  VectorXd newmEffcts = mEffcts;
  newmEffcts(proposal.mTile) = proposal.newmEffct;
  calc_B(nowmColors,newmEffcts,mrateMu,proposal.newB);
  return (EEMS_wishpdfln(proposal.newB,nowq,s2loc,proposal.newtrDinvQxD));
}
double EEMS::eval_proposal_mrateMu(Proposal &proposal) const {
  calc_B(nowmColors,mEffcts,proposal.newmrateMu,proposal.newB);
  return (EEMS_wishpdfln(proposal.newB,nowq,s2loc,proposal.newtrDinvQxD));
}
double EEMS::eval_proposal_mSeeds(Proposal &proposal) const {
  MatrixXd newmSeeds = mSeeds;
  newmSeeds(proposal.mTile,0) = proposal.newmSeedx;
  newmSeeds(proposal.mTile,1) = proposal.newmSeedy;
  graph.pdist2(newmSeeds,proposal.newmColors);
  calc_B(proposal.newmColors,mEffcts,mrateMu,proposal.newB);
  return (EEMS_wishpdfln(proposal.newB,nowq,s2loc,proposal.newtrDinvQxD));
}
double EEMS::eval_birthdeath_mVoronoi(Proposal &proposal) const {
  graph.pdist2(proposal.newmSeeds,proposal.newmColors);
  calc_B(proposal.newmColors,proposal.newmEffcts,mrateMu,proposal.newB);
  return (EEMS_wishpdfln(proposal.newB,nowq,s2loc,proposal.newtrDinvQxD));
}
///////////////////////////////////////////
// Proposals:
void EEMS::update_s2loc( ) {
  for ( int i = 0 ; i < p ; i++ ) {
    nowpi += (params.s2locShape/2.0+1.0)*log(s2loc(i)) + (params.s2locScale/2.0)/s2loc(i);
    nowll += 0.5 * nowtrDinvQxD(i) /s2loc(i) + 0.5 * nmin1(i) * log(s2loc(i));
    s2loc(i) = draw.rinvgam( 0.5 * ( params.s2locShape + nmin1(i) ),
			     0.5 * ( params.s2locScale + nowtrDinvQxD(i) ) );
    nowpi -= (params.s2locShape/2.0+1.0)*log(s2loc(i)) + (params.s2locScale/2.0)/s2loc(i);
    nowll -= 0.5 * nowtrDinvQxD(i) /s2loc(i) + 0.5 * nmin1(i) * log(s2loc(i));
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
      nowB = proposal.newB;
      break;
    case 4:
      mrateMu = proposal.newmrateMu;
      nowB = proposal.newB;
      break;
    case 5:
      mTile = proposal.mTile;
      mSeeds(mTile,0) = proposal.newmSeedx;
      mSeeds(mTile,1) = proposal.newmSeedy;
      nowB = proposal.newB;
      nowmColors = proposal.newmColors;
      break;
    default:
      mSeeds = proposal.newmSeeds;
      mtiles = proposal.newmtiles;
      mEffcts = proposal.newmEffcts;
      nowB = proposal.newB;
      nowmColors = proposal.newmColors;
    }
    nowpi = proposal.newpi;
    nowll = proposal.newll;
    nowtrDinvQxD = proposal.newtrDinvQxD;
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
  mcmcthetas.row(iter) = s2loc;
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
  JtDhatJ_allSites += nowB + 0.5 * nowq.replicate(1,o) + 0.5 * nowq.transpose().replicate(o,1);
  return true;
}
bool EEMS::output_results(const MCMC &mcmc) const {
  ofstream out; bool done;
  MatrixXd oDemes = MatrixXd::Zero(o,3);
  oDemes << graph.get_the_obsrv_demes(),c_allSites;
  out.open((params.mcmcpath + "/rdistoDemes.txt").c_str(),ofstream::out);
  if (!out.is_open()) { return false; }
  out << oDemes << endl;
  out.close( );
  out.open((params.mcmcpath + "/rdistJtDobsJ.txt").c_str(),ofstream::out);
  if (!out.is_open()) { return false; }
  MatrixXd Pairs_allSites = c_allSites*c_allSites.transpose();
  Pairs_allSites -= c_allSites.asDiagonal();
  out << JtDobsJ_allSites.cwiseQuotient(Pairs_allSites) << endl;
  out.close( );
  out.open((params.mcmcpath + "/rdistJtDhatJ.txt").c_str(),ofstream::out);
  if (!out.is_open()) { return false; }
  out << JtDhatJ_allSites/niters << endl;
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
  done = dlmcell(params.mcmcpath + "/mcmcmrates.txt",mcmcmtiles,mcmcmRates); if (!done) { return false; }
  done = dlmcell(params.mcmcpath + "/mcmcxcoord.txt",mcmcmtiles,mcmcxCoord); if (!done) { return false; }
  done = dlmcell(params.mcmcpath + "/mcmcycoord.txt",mcmcmtiles,mcmcyCoord); if (!done) { return false; }
  done = dlmcell(params.mcmcpath + "/mcmcqrates.txt",mcmcqtiles,mcmcqRates); if (!done) { return false; }
  done = dlmcell(params.mcmcpath + "/mcmcwcoord.txt",mcmcqtiles,mcmcwCoord); if (!done) { return false; }
  done = dlmcell(params.mcmcpath + "/mcmczcoord.txt",mcmcqtiles,mcmczCoord); if (!done) { return false; }
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
       << "and effective degrees of freedom = " << p << endl
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
string EEMS::datapath( ) const { return params.datapath; }
string EEMS::mcmcpath( ) const { return params.mcmcpath; }
