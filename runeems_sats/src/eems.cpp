
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
    Bconst = 1.0;  qconst = 2.0;
  } else {
    Bconst = 0.25; qconst = 1.0;
  }
}
EEMS::~EEMS( ) { }
///////////////////////////////////////////
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
    int nimin1 = ni-1;
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
    MatrixXd Li = MatrixXd::Constant(nimin1,ni,-1.0);
    Li.topRightCorner(nimin1,nimin1).setIdentity();
    Z.push_back(zi);
    O.push_back(oi);
    L.push_back(Li);
    // Similarity matrix (with rank 1)
    MatrixXd Si = zi*zi.transpose();
    // Dissimilarity matrix
    MatrixXd Di = - 2.0 * Si;
    Di += Si.diagonal().replicate(1,ni);
    Di += Si.diagonal().transpose().replicate(ni,1);
    ll_partdf += logdet(Li*Li.transpose()) + nimin1 * (log_2 + log_pi)
      + nimin1 * pseudologdet( - Li*Di*Li.transpose(),1);
    logn(i) = log(ni);
    nmin1(i) = nimin1;
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
void EEMS::initialize_state( ) {
  cerr << "[EEMS::initialize_state]" << endl;
  sigma2 = VectorXd::Zero(p);
  for (int i = 0 ; i < p ; i++ ) {
    sigma2(i) = draw.rinvgam(3.0,1.0);
  }
  // Initialize the two Voronoi tessellations
  qtiles = draw.rnegbin(2*o,0.5); // o is the number of observed demes
  mtiles = draw.rnegbin(2*o,0.5);
  cerr << "  EEMS starts with " << qtiles << " qtiles and " << mtiles << " mtiles" << endl;
  // Draw the Voronoi centers Coord uniformly within the habitat
  qSeeds = MatrixXd::Zero(qtiles,2); randpoint_in_habitat(qSeeds);
  mSeeds = MatrixXd::Zero(mtiles,2); randpoint_in_habitat(mSeeds);
  mrateS2 = draw.rinvgam(0.5,0.5);
  qrateS2 = draw.rinvgam(0.5,0.5);
  // Assign migration rates to the Voronoi tiles
  mrateMu = params.mrateMuHalfInterval*(2.0*draw.runif() - 1.0);
  // Assign rates to the Voronoi tiles
  qEffcts = VectorXd::Zero(qtiles); rnorm_effects(params.qEffctHalfInterval,qrateS2,qEffcts);
  mEffcts = VectorXd::Zero(mtiles); rnorm_effects(params.mEffctHalfInterval,mrateS2,mEffcts);
  cerr << "[EEMS::initialize_state] Done." << endl << endl;
}
void EEMS::load_final_state( ) {
  cerr << "[EEMS::load_final_state]" << endl;  
  MatrixXd tempi; bool error = false;
  tempi = readMatrixXd(params.prevpath + "/lastqtiles.txt");
  if ((tempi.rows()!=1) || (tempi.cols()!=1)) { error = true; }
  qtiles = tempi(0,0);
  tempi = readMatrixXd(params.prevpath + "/lastmtiles.txt");
  if ((tempi.rows()!=1) || (tempi.cols()!=1)) { error = true; }
  mtiles = tempi(0,0);
  cerr << "  EEMS starts with " << qtiles << " qtiles and " << mtiles << " mtiles" << endl;
  tempi = readMatrixXd(params.prevpath + "/lastthetas.txt");
  if ((tempi.rows()!=p) || (tempi.cols()!=1)) { error = true; }
  sigma2 = tempi.col(0);
  tempi = readMatrixXd(params.prevpath + "/lastqhyper.txt");
  if ((tempi.rows()!=1) || (tempi.cols()!=1)) { error = true; }
  qrateS2 = tempi(0,0);
  tempi = readMatrixXd(params.prevpath + "/lastmhyper.txt");
  if ((tempi.rows()!=1) || (tempi.cols()!=2)) { error = true; }
  mrateMu = tempi(0,0);
  mrateS2 = tempi(0,1);
  tempi = readMatrixXd(params.prevpath + "/lastqeffct.txt");
  if ((tempi.rows()!=qtiles) || (tempi.cols()!=1)) { error = true; }
  qEffcts = tempi.col(0);
  tempi = readMatrixXd(params.prevpath + "/lastmeffct.txt");
  if ((tempi.rows()!=mtiles) || (tempi.cols()!=1)) { error = true; }
  mEffcts = tempi.col(0);
  qSeeds = readMatrixXd(params.prevpath + "/lastqseeds.txt");
  if ((qSeeds.rows()!=qtiles) || (qSeeds.cols()!=2)) { error = true; }
  mSeeds = readMatrixXd(params.prevpath + "/lastmseeds.txt");
  if ((mSeeds.rows()!=mtiles) || (mSeeds.cols()!=2)) { error = true; }
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
    move = M_MEAN_RATE_UPDATE;
  }
  return(move);
}
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
  if (!inrange) { return (-Inf); }
  nowpi = 
    + lgamma(params.negBiSize+mtiles) - lgamma(mtiles+1.0) + mtiles*log(params.negBiProb)
    + lgamma(params.negBiSize+qtiles) - lgamma(qtiles+1.0) + qtiles*log(params.negBiProb)
    - (params.mrateShape_2+1.0)*log(mrateS2) - params.mrateScale_2/mrateS2
    - (params.qrateShape_2+1.0)*log(qrateS2) - params.qrateScale_2/qrateS2
    - (mtiles/2.0)*log(mrateS2) - mEffcts.squaredNorm()/(2.0*mrateS2)
    - (qtiles/2.0)*log(qrateS2) - qEffcts.squaredNorm()/(2.0*qrateS2)
    - (params.sigmaShape_2+1.0)*sigma2.array().log().sum()
    - params.sigmaScale_2*pow(sigma2.array(),-1.0).sum();
  return (nowpi);
}
double EEMS::eval_likelihood( ) {
  graph.index_closest_to_deme(mSeeds,nowmColors);
  graph.index_closest_to_deme(qSeeds,nowqColors);
  calc_q(nowqColors,qEffcts,nowq);
  calc_B(nowmColors,mEffcts,mrateMu,nowB);
  VectorXd trDinvQxD = VectorXd::Zero(p);
  nowll = EEMS_wishpdfln(nowB,nowq,sigma2,trDinvQxD);
  nowtrDinvQxD = trDinvQxD;
  return (nowll);
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
double EEMS::EEMS_wishpdfln(const MatrixXd &B, const VectorXd &q, const VectorXd &sigma2, VectorXd &trDinvQxD) const {
  VectorXd ldetDinvQ = VectorXd::Zero(p);
  if (trDinvQxD.size() != p) { trDinvQxD.resize(p); }
  for ( int i = 0 ; i < p ; i++ ) {
    VectorXd qi = slice(q,ovec[i]);
    VectorXd qiinv = pow(qi.array(),-1.0);
    MatrixXd T = slice(B,ovec[i],ovec[i]);
    T *= cvec[i].asDiagonal(); T -= qi.asDiagonal(); // Now T = B*C - W
    PartialPivLU<MatrixXd> lu(T);
    T.noalias() = T*qiinv.asDiagonal()*cinv[i].asDiagonal();
    T += cinv[i].asDiagonal();                       // Now T = B*inv(W)
    MatrixXd X = lu.solve(T);
    VectorXd Xc_qinv = X*cvec[i] - qiinv;
    double oDinvo = cvec[i].dot(Xc_qinv);
    double oDiDDi = Xc_qinv.transpose()*JtDobsJ[i]*Xc_qinv;
    ldetDinvQ(i) = logn(i) - log(abs(oDinvo))
      + cmin1[i].dot(qiinv.array().log().matrix())
      - lu.matrixLU().diagonal().array().abs().log().sum();
    trDinvQxD(i) = trace_AxB(X,JtDobsJ[i]) - oDiDDi/oDinvo;
  }
  return ( 0.5 * ( ldetDinvQ.sum() - (trDinvQxD.array() / sigma2.array()).sum() -
		   (nmin1.array() * sigma2.array().log()).sum() - ll_partdf ) );
}
double EEMS::eval_proposal_qEffcts(Proposal &proposal) const {
  VectorXd newqEffcts = qEffcts;
  newqEffcts(proposal.qTile) = proposal.newqEffct;
  calc_q(nowqColors,newqEffcts,proposal.newq);
  return (EEMS_wishpdfln(nowB,proposal.newq,sigma2,proposal.newtrDinvQxD));
}
double EEMS::eval_proposal_qSeeds(Proposal &proposal) const {
  MatrixXd newqSeeds = qSeeds;
  newqSeeds(proposal.qTile,0) = proposal.newqSeedx;
  newqSeeds(proposal.qTile,1) = proposal.newqSeedy;
  graph.index_closest_to_deme(newqSeeds,proposal.newqColors);
  calc_q(proposal.newqColors,qEffcts,proposal.newq);
  return (EEMS_wishpdfln(nowB,proposal.newq,sigma2,proposal.newtrDinvQxD));
}
double EEMS::eval_birthdeath_qVoronoi(Proposal &proposal) const {
  graph.index_closest_to_deme(proposal.newqSeeds,proposal.newqColors);
  calc_q(proposal.newqColors,proposal.newqEffcts,proposal.newq);
  return (EEMS_wishpdfln(nowB,proposal.newq,sigma2,proposal.newtrDinvQxD));
}
double EEMS::eval_proposal_mEffcts(Proposal &proposal) const {
  VectorXd newmEffcts = mEffcts;
  newmEffcts(proposal.mTile) = proposal.newmEffct;
  calc_B(nowmColors,newmEffcts,mrateMu,proposal.newB);
  return (EEMS_wishpdfln(proposal.newB,nowq,sigma2,proposal.newtrDinvQxD));
}
double EEMS::eval_proposal_mrateMu(Proposal &proposal) const {
  calc_B(nowmColors,mEffcts,proposal.newmrateMu,proposal.newB);
  return (EEMS_wishpdfln(proposal.newB,nowq,sigma2,proposal.newtrDinvQxD));
}
double EEMS::eval_proposal_mSeeds(Proposal &proposal) const {
  MatrixXd newmSeeds = mSeeds;
  newmSeeds(proposal.mTile,0) = proposal.newmSeedx;
  newmSeeds(proposal.mTile,1) = proposal.newmSeedy;
  graph.index_closest_to_deme(newmSeeds,proposal.newmColors);
  calc_B(proposal.newmColors,mEffcts,mrateMu,proposal.newB);
  return (EEMS_wishpdfln(proposal.newB,nowq,sigma2,proposal.newtrDinvQxD));
}
double EEMS::eval_birthdeath_mVoronoi(Proposal &proposal) const {
  graph.index_closest_to_deme(proposal.newmSeeds,proposal.newmColors);
  calc_B(proposal.newmColors,proposal.newmEffcts,mrateMu,proposal.newB);
  return (EEMS_wishpdfln(proposal.newB,nowq,sigma2,proposal.newtrDinvQxD));
}
///////////////////////////////////////////
void EEMS::update_sigma2( ) {
  for ( int i = 0 ; i < p ; i++ ) {
    nowpi += (params.sigmaShape_2+1.0)*log(sigma2(i)) + params.sigmaScale_2/sigma2(i);
    nowll += 0.5 * nowtrDinvQxD(i) /sigma2(i) + 0.5 * nmin1(i) * log(sigma2(i));
    sigma2(i) = draw.rinvgam( params.sigmaShape_2 + 0.5 * nmin1(i),
			      params.sigmaScale_2 + 0.5 * nowtrDinvQxD(i) );
    nowpi -= (params.sigmaShape_2+1.0)*log(sigma2(i)) + params.sigmaScale_2/sigma2(i);
    nowll -= 0.5 * nowtrDinvQxD(i) /sigma2(i) + 0.5 * nmin1(i) * log(sigma2(i));
  }
}
void EEMS::propose_qEffcts(Proposal &proposal) {
  // Choose a tile at random to update
  proposal.qTile = draw.runif_int(0,qtiles-1);
  double curqEffct = qEffcts(proposal.qTile);
  double newqEffct = draw.rnorm(curqEffct,params.qEffctProposalS2);
  proposal.type = 0;
  proposal.newqEffct = newqEffct;
  // The prior distribution on the tile effects is truncated normal
  // So first check whether the proposed value is in range
  // Then update the prior and evaluate the new likelihood
  //  + (curqEffct*curqEffct) / (2.0*qrateS2) : old prior component associated with this tile
  //  - (newqEffct*newqEffct) / (2.0*qrateS2) : new prior component associated with this tile
  if ( abs(newqEffct) < params.qEffctHalfInterval ) {
    proposal.newpi = nowpi - (newqEffct*newqEffct - curqEffct*curqEffct) / (2.0*qrateS2);
    proposal.newll = eval_proposal_qEffcts(proposal);
  } else {
    proposal.newpi = -Inf;
    proposal.newll = -Inf;
  }
}
void EEMS::propose_mEffcts(Proposal &proposal) {
  proposal.mTile = draw.runif_int(0,mtiles-1);
  double curmEffct = mEffcts(proposal.mTile);
  double newmEffct = draw.rnorm(curmEffct,params.mEffctProposalS2);
  proposal.type = 3;
  proposal.newmEffct = newmEffct;
  if ( abs(newmEffct) < params.mEffctHalfInterval ) {
    proposal.newpi = nowpi - (newmEffct*newmEffct - curmEffct*curmEffct) / (2.0*mrateS2);
    proposal.newll = eval_proposal_mEffcts(proposal);
  } else {
    proposal.newpi = -Inf;
    proposal.newll = -Inf;
  }
}
void EEMS::propose_mrateMu(Proposal &proposal) {
  // Make a random-walk Metropolis-Hastings proposal
  double newmrateMu = draw.rnorm(mrateMu,params.mrateMuProposalS2);
  proposal.type = 4;
  proposal.newmrateMu = newmrateMu;
  // If the proposed value is in range, the prior probability does not change
  // as the prior distribution on mrateMu is uniform
  // Otherwise, setting the prior and the likelihood to -Inf forces a rejection
  if ( abs(newmrateMu) < params.mrateMuHalfInterval ) {
    proposal.newpi = nowpi;
    proposal.newll = eval_proposal_mrateMu(proposal);
  } else {
    proposal.newpi = -Inf;
    proposal.newll = -Inf;
  }
}
void EEMS::move_qVoronoi(Proposal &proposal) {
  // Choose a tile at random to move
  proposal.qTile = draw.runif_int(0,qtiles-1);
  // A move is a jitter about the current (x,y) location of the tile center (seed)
  double newqSeedx = draw.rnorm(qSeeds(proposal.qTile,0),params.qSeedsProposalS2x);
  double newqSeedy = draw.rnorm(qSeeds(proposal.qTile,1),params.qSeedsProposalS2y);
  proposal.type = 1;
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
void EEMS::move_mVoronoi(Proposal &proposal) {
  proposal.mTile = draw.runif_int(0,mtiles-1);
  double newmSeedx = draw.rnorm(mSeeds(proposal.mTile,0),params.mSeedsProposalS2x);
  double newmSeedy = draw.rnorm(mSeeds(proposal.mTile,1),params.mSeedsProposalS2y);
  proposal.type = 5;
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
void EEMS::birthdeath_qVoronoi(Proposal &proposal) {
  int newqtiles = qtiles,r;
  double u = draw.runif();
  double pBirth = 0.5;
  double pDeath = 0.5;
  boost::math::normal pnorm(0.0,sqrt(qrateS2));
  proposal.newqEffcts = qEffcts;
  proposal.newqSeeds = qSeeds;
  // If there is exactly one tile, rule out a death proposal
  if ((qtiles==1) || (u<0.5)) { // Propose birth
    if (qtiles==1) { pBirth = 1.0; }
    newqtiles++;
    MatrixXd newqSeed = MatrixXd::Zero(1,2);
    randpoint_in_habitat(newqSeed);
    pairwise_distance(qSeeds,newqSeed).col(0).minCoeff(&r);
    // The new tile is assigned a rate by perturbing the current rate at the new seed    
    double nowqEffct = qEffcts(r);
    double newqEffct = draw.rtrnorm(nowqEffct,params.qEffctProposalS2,params.qEffctHalfInterval);
    insertRow(proposal.newqSeeds,newqSeed.row(0));
    insertElem(proposal.newqEffcts,newqEffct);
    // Compute log(proposal ratio) and log(prior ratio)
    proposal.ratioln = log(pDeath/pBirth)
      - dtrnormln(newqEffct,nowqEffct,params.qEffctProposalS2,params.qEffctHalfInterval)
      - log(cdf(pnorm,params.qEffctHalfInterval) - cdf(pnorm,-params.qEffctHalfInterval));
    proposal.newpi = nowpi + log((qtiles+params.negBiSize)/(newqtiles/params.negBiProb))
      - 0.5 * log(qrateS2) - 0.5 * newqEffct*newqEffct/qrateS2;
  } else {                      // Propose death
    if (qtiles==2) { pBirth = 1.0; }
    newqtiles--;
    int qtileToRemove = draw.runif_int(0,newqtiles);
    MatrixXd oldqSeed = qSeeds.row(qtileToRemove);
    removeRow(proposal.newqSeeds,qtileToRemove);
    removeElem(proposal.newqEffcts,qtileToRemove);
    pairwise_distance(proposal.newqSeeds,oldqSeed).col(0).minCoeff(&r);
    double nowqEffct = proposal.newqEffcts(r);
    double oldqEffct = qEffcts(qtileToRemove);
    // Compute log(prior ratio) and log(proposal ratio)
    proposal.ratioln = log(pBirth/pDeath)
      + dtrnormln(oldqEffct,nowqEffct,params.qEffctProposalS2,params.qEffctHalfInterval)
      + log(cdf(pnorm,params.qEffctHalfInterval) - cdf(pnorm,-params.qEffctHalfInterval));
    proposal.newpi = nowpi + log((qtiles/params.negBiProb)/(newqtiles+params.negBiSize))
      + 0.5 * log(qrateS2) + 0.5 * oldqEffct*oldqEffct/qrateS2;
  }
  proposal.type = 2;
  proposal.newqtiles = newqtiles;
  proposal.newll = eval_birthdeath_qVoronoi(proposal);
}
void EEMS::birthdeath_mVoronoi(Proposal &proposal) {
  int newmtiles = mtiles,r;
  double u = draw.runif();
  double pBirth = 0.5;
  double pDeath = 0.5;
  boost::math::normal pnorm(0.0,sqrt(mrateS2));
  proposal.newmEffcts = mEffcts;
  proposal.newmSeeds = mSeeds;
  if ((mtiles==1) || (u<0.5)) { // Propose birth
    if (mtiles==1) { pBirth = 1.0; }
    newmtiles++;
    MatrixXd newmSeed = MatrixXd::Zero(1,2);
    randpoint_in_habitat(newmSeed);
    pairwise_distance(mSeeds,newmSeed).col(0).minCoeff(&r);
    double nowmEffct = mEffcts(r);
    double newmEffct = draw.rtrnorm(nowmEffct,params.mEffctProposalS2,params.mEffctHalfInterval);
    insertRow(proposal.newmSeeds,newmSeed.row(0));
    insertElem(proposal.newmEffcts,newmEffct);
    proposal.ratioln = log(pDeath/pBirth)
      - dtrnormln(newmEffct,nowmEffct,params.mEffctProposalS2,params.mEffctHalfInterval)
      - log(cdf(pnorm,params.mEffctHalfInterval) - cdf(pnorm,-params.mEffctHalfInterval));
    proposal.newpi = nowpi + log((mtiles+params.negBiSize)/(newmtiles/params.negBiProb))
      - 0.5 * log(mrateS2) - 0.5 * newmEffct*newmEffct/mrateS2;
  } else {                      // Propose death
    if (mtiles==2) { pBirth = 1.0; }
    newmtiles--;
    int mtileToRemove = draw.runif_int(0,newmtiles);
    MatrixXd oldmSeed = mSeeds.row(mtileToRemove);
    removeRow(proposal.newmSeeds,mtileToRemove);
    removeElem(proposal.newmEffcts,mtileToRemove);
    pairwise_distance(proposal.newmSeeds,oldmSeed).col(0).minCoeff(&r);
    double nowmEffct = proposal.newmEffcts(r);
    double oldmEffct = mEffcts(mtileToRemove);
    proposal.ratioln = log(pBirth/pDeath)
      + dtrnormln(oldmEffct,nowmEffct,params.mEffctProposalS2,params.mEffctHalfInterval)      
      + log(cdf(pnorm,params.mEffctHalfInterval) - cdf(pnorm,-params.mEffctHalfInterval));
    proposal.newpi = nowpi + log((mtiles/params.negBiProb)/(newmtiles+params.negBiSize))
      + 0.5 * log(mrateS2) + 0.5 * oldmEffct*oldmEffct/mrateS2;
  }
  proposal.type = 6;
  proposal.newmtiles = newmtiles;
  proposal.newll = eval_birthdeath_mVoronoi(proposal);
}
void EEMS::update_hyperparams( ) {
  double SSq = qEffcts.squaredNorm();
  double SSm = mEffcts.squaredNorm();
  nowpi += (params.mrateShape_2+1.0)*log(mrateS2) + params.mrateScale_2/mrateS2
    +      (params.qrateShape_2+1.0)*log(qrateS2) + params.qrateScale_2/qrateS2
    + (mtiles/2.0)*log(mrateS2) + SSm/(2.0*mrateS2)
    + (qtiles/2.0)*log(qrateS2) + SSq/(2.0*qrateS2);
  qrateS2 = draw.rinvgam(params.qrateShape_2 + 0.5 * qtiles, params.qrateScale_2 + 0.5 * SSq);
  mrateS2 = draw.rinvgam(params.mrateShape_2 + 0.5 * mtiles, params.mrateScale_2 + 0.5 * SSm);
  nowpi -= (params.mrateShape_2+1.0)*log(mrateS2) + params.mrateScale_2/mrateS2
    +      (params.qrateShape_2+1.0)*log(qrateS2) + params.qrateScale_2/qrateS2
    + (mtiles/2.0)*log(mrateS2) + SSm/(2.0*mrateS2)
    + (qtiles/2.0)*log(qrateS2) + SSq/(2.0*qrateS2);
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
  if (proposal.type==2 || proposal.type==6) {
    ratioln += proposal.ratioln;
  }
  if ( log(u) < min(0.0,ratioln) ) {
    switch (proposal.type) {
    case 0:
      qEffcts(proposal.qTile) = proposal.newqEffct;
      nowq = proposal.newq;
      break;
    case 1:
      qSeeds(proposal.qTile,0) = proposal.newqSeedx;
      qSeeds(proposal.qTile,1) = proposal.newqSeedy;
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
      mEffcts(proposal.mTile) = proposal.newmEffct;
      nowB = proposal.newB;
      break;
    case 4:
      mrateMu = proposal.newmrateMu;
      nowB = proposal.newB;
      break;
    case 5:
      mSeeds(proposal.mTile,0) = proposal.newmSeedx;
      mSeeds(proposal.mTile,1) = proposal.newmSeedy;
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
       << "         number of qVoronoi tiles = " << qtiles << endl
       << "         number of mVoronoi tiles = " << mtiles << endl
       << "          Log prior = " << nowpi << endl
       << "          Log llike = " << nowll << endl;
}
void EEMS::save_iteration(const MCMC &mcmc) {
  int iter = mcmc.to_save_iteration( );
  mcmcthetas.row(iter) = sigma2;
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
}
bool EEMS::output_current_state( ) const {
  ofstream out; bool error = false;
  out.open((params.mcmcpath + "/lastqtiles.txt").c_str(),ofstream::out);
  if (!out.is_open()) { error = true; return(error); }
  out << qtiles << endl;
  out.close( );
  out.open((params.mcmcpath + "/lastmtiles.txt").c_str(),ofstream::out);
  if (!out.is_open()) { error = true; return(error); }
  out << mtiles << endl;
  out.close( );
  out.open((params.mcmcpath + "/lastthetas.txt").c_str(),ofstream::out);
  if (!out.is_open()) { error = true; return(error); }
  out << fixed << setprecision(6) << sigma2 << endl;
  out.close( );
  out.open((params.mcmcpath + "/lastqhyper.txt").c_str(),ofstream::out);
  if (!out.is_open()) { error = true; return(error); }
  out << fixed << setprecision(6) << qrateS2 << endl;
  out.close( );
  out.open((params.mcmcpath + "/lastmhyper.txt").c_str(),ofstream::out);
  if (!out.is_open()) { error = true; return(error); }
  out << fixed << setprecision(6) << mrateMu << " " << mrateS2 << endl;
  out.close( );
  out.open((params.mcmcpath + "/lastpilogl.txt").c_str(),ofstream::out);
  if (!out.is_open()) { error = true; return(error); }
  out << fixed << setprecision(6) << nowpi << " " << nowll << endl;
  out.close( );
  out.open((params.mcmcpath + "/lastmeffct.txt").c_str(),ofstream::out);
  if (!out.is_open()) { error = true; return(error); }
  out << fixed << setprecision(6) << mEffcts << endl;
  out.close( );
  out.open((params.mcmcpath + "/lastmseeds.txt").c_str(),ofstream::out);
  if (!out.is_open()) { error = true; return(error); }
  out << fixed << setprecision(6) << mSeeds << endl;
  out.close( );
  out.open((params.mcmcpath + "/lastqeffct.txt").c_str(),ofstream::out);
  if (!out.is_open()) { error = true; return(error); }
  out << fixed << setprecision(6) << qEffcts << endl;
  out.close( );
  out.open((params.mcmcpath + "/lastqseeds.txt").c_str(),ofstream::out);
  if (!out.is_open()) { error = true; return(error); }
  out << fixed << setprecision(6) << qSeeds << endl;
  out.close( );
  return(error);
}
bool EEMS::output_results(const MCMC &mcmc) const {
  ofstream out; bool error = false;
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
  int niters = mcmc.num_iters_to_save();
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
  double pi0 = test_prior( );
  double ll0 = test_likelihood( );
  if ((abs(nowpi-pi0)/abs(pi0)>1e-12)||
      (abs(nowll-ll0)/abs(ll0)>1e-12)) {
    cerr << "[EEMS::testing]   |ll0-ll|/|ll0| = " << abs(nowll - ll0)/abs(ll0) << endl;
    cerr << "[EEMS::testing]   |pi0-pi|/|pi0| = " << abs(nowpi - pi0)/abs(pi0) << endl;
    exit(1);
  }
}
double EEMS::test_prior( ) const {
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
  if (!inrange) { return (-Inf); }
  double logpi =
    + dnegbinln(mtiles,params.negBiSize,params.negBiProb)
    + dnegbinln(qtiles,params.negBiSize,params.negBiProb)
    + dinvgamln(mrateS2,params.mrateShape_2,params.mrateScale_2)
    + dinvgamln(qrateS2,params.qrateShape_2,params.qrateScale_2)
    + dmvnormln(mEffcts,VectorXd::Zero(mtiles),mrateS2*MatrixXd::Identity(mtiles,mtiles))
    + dmvnormln(qEffcts,VectorXd::Zero(qtiles),qrateS2*MatrixXd::Identity(qtiles,qtiles));
  for ( int i = 0 ; i < p ; i++ ) {
    logpi += dinvgamln(sigma2(i),params.sigmaShape_2,params.sigmaScale_2);
  }
  return (logpi);
}
double EEMS::test_likelihood( ) const {
  VectorXi mColors, qColors;
  graph.index_closest_to_deme(mSeeds,mColors);
  graph.index_closest_to_deme(qSeeds,qColors);
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
			     -L[i]*slice(Delta,O[i],O[i])*L[i].transpose()*sigma2(i),1);
  }
  return (logll);
}
