#pragma once

#include "util.hpp"
#include "mcmc.hpp"
#include "draw.hpp"
#include "graph.hpp"
#include "habitat.hpp"

#include <boost/filesystem.hpp>

#ifndef EEMS_H
#define EEMS_H

/*
  An updated set of parameter values
  The type of move is necessary in order to know which parameters have a new proposed value;
  the rest of the parameters won't be set to their current values (to avoid unnecessary copying)
  For example, if move = M_VORONOI_BIRTH_DEATH,
               then newmtiles, newmSeeds, nowmColors, newmEffcts (and of course, newpi, newll, newratioln) would have changed
	       if move = M_VORONOI_POINT_MOVE,
	       then newmSeeds, nowmColors (and of course, newpi, newll) would have changed
  The ratioln is the proposal ratio for birth/death proposal.
  For the usual Metropolis-Hastings updates, the acceptance probability is
    alpha = (prior ratio) * (likelihood ratio)
  For the birth/deatch RJ-MCMC updates, the acceptance probability is
    alpha = (proposal ratio) * (prior ratio) * (likelihood ratio)
  See Green, "Reversible jump Markov chain Monte Carlo computation and Bayesian model determination"
*/
struct Proposal {
  MoveType move; // the type of proposal/update
  int newqtiles; // number of m and q tiles, respectively
  int newmtiles;
  double newpi; // log prior
  double newll; // log likelihood
  double newsigma2; // variance scale
  double newratioln; // RJ-MCMC proposal ratio, on the log scale
  double newmrateMu; // overall (mean) migration rate, on the log10 scale; qrateMu is assumed to be 0
  VectorXd newtriDeltaQD; // precomputed components of the Wishart log likelihood, one for each microsatellite locus
  VectorXd newqEffcts; // the diversity rate of each q tile, on the log10 scale
  VectorXd newmEffcts; // the migration rate of each m tile, on the log10 scale and relative to the ovarall mrateMu
  MatrixXd newqSeeds;  // the location of each q tile within the habitat
  MatrixXd newmSeeds;  // the location of each m tile within the habitat
  VectorXd newW; // the within demes component, one entry for each deme in the habitat
  MatrixXd newB; // the between demes component, one entry for each pair of demes in the habitat
  VectorXi newqColors; // mapping that indicates which q tiles each vertex/deme falls into
  VectorXi newmColors; // mapping that indicates which m tiles each vertex/deme falls into
};

class EEMS {
public:

  EEMS(const Params &params);
  ~EEMS( );

  void initialize_state( );
  void load_final_state( );
  bool start_eems(const MCMC &mcmc);
  double eval_prior( );
  double eval_likelihood( );
  double test_prior(const MatrixXd &mSeeds, const VectorXd &mEffcts, const double mrateMu,
		    const MatrixXd &qSeeds, const VectorXd &qEffcts,
		    const VectorXd &sigma2, const double mrateS2, const double qrateS2) const;
  double test_likelihood(const MatrixXd &mSeeds, const VectorXd &mEffcts, const double mrateMu,
			 const MatrixXd &qSeeds, const VectorXd &qEffcts,
			 const VectorXd &sigma2) const;
  void calc_within(const VectorXi &qColors, const VectorXd &qEffcts, VectorXd &W) const;
  void calc_between(const VectorXi &mColors, const VectorXd &mEffcts, const double mrateMu, MatrixXd &B) const;
  MoveType choose_move_type( );
  // These functions change the within demes component:
  double eval_proposal_rate_one_qtile(Proposal &proposal) const;
  double eval_proposal_move_one_qtile(Proposal &proposal) const;
  double eval_birthdeath_qVoronoi(Proposal &proposal) const;
  // These functions change the between demes component:
  double eval_proposal_rate_one_mtile(Proposal &proposal) const;
  double eval_proposal_overall_mrate(Proposal &proposal) const;
  double eval_proposal_move_one_mtile(Proposal &proposal) const;
  double eval_birthdeath_mVoronoi(Proposal &proposal) const;

  // Gibbs updates:
  void update_sigma2( );
  void update_hyperparams( );
  // Random-walk Metropolis-Hastings proposals:
  void propose_df(Proposal &proposal,const MCMC &mcmc);
  void propose_rate_one_qtile(Proposal &proposal);
  void propose_rate_one_mtile(Proposal &proposal);
  void propose_overall_mrate(Proposal &proposal);
  void propose_move_one_qtile(Proposal &proposal);
  void propose_move_one_mtile(Proposal &proposal);
  void propose_birthdeath_qVoronoi(Proposal &proposal);
  void propose_birthdeath_mVoronoi(Proposal &proposal);
  bool accept_proposal(Proposal &proposal);
  ///////////////////////////////////////////

  void print_iteration(const MCMC &mcmc) const;
  void save_iteration(const MCMC &mcmc);
  void output_results(const MCMC &mcmc) const;
  void output_current_state() const;
  void check_ll_computation() const;
  string datapath() const;
  string mcmcpath() const;
  string prevpath() const;
  string gridpath() const;

private:

  Draw draw; // Random number generator
  Graph graph;
  Params params;
  Habitat habitat;

  // Diffs:
  int o; // observed demes
  int d; // all demes
  int n; // nIndiv
  int p; // nSites
  MatrixXd J;
  vector<MatrixXd> Diffs;
  vector<MatrixXd> L;
  vector<MatrixXd> Z;
  vector<MatrixXi> O;
  vector<VectorXi> ovec;
  // cvec, cinv, cmin1 are used to compute EEMS_wishpdfln
  vector<VectorXd> cvec; // c is the vector of counts
  vector<VectorXd> cinv;  // cinv is the vector of inverse counts
  vector<VectorXd> cmin1;  // cmin1 is the vector of counts - 1
  vector<MatrixXd> JtDobsJ;
  VectorXd c_allSites;
  MatrixXd JtDobsJ_allSites;
  MatrixXd JtDhatJ_allSites;
  double ldLLt; // logdet(L*L')
  double ldDiQ;  // logdet(inv(Diffs)*Q)
  double ldLDLt;  // logdet(-L*Diffs*L')
  VectorXd logn;
  VectorXd nmin1;
  void initialize_diffs();
  void randpoint_in_habitat(MatrixXd &Seeds);
  void rnorm_effects(const double HalfInterval, const double rateS2, VectorXd &Effcts);
  
  // The current set of parameter values:
  int nowmtiles, nowqtiles; // number of m and q tiles, respectively
  MatrixXd nowmSeeds; VectorXd nowmEffcts; double nowmrateMu; // parameters to describe the m Voronoi tessellation
  MatrixXd nowqSeeds; VectorXd nowqEffcts;                    // parameters to describe the q Voronoi tessellation
  VectorXd nowsigma2; // a vector of scale variances, one for each microsatellite locus
  double nowqrateS2, nowmrateS2; // two hyperparameters -- the variance of nowqEffcts and nowmEffcts, respectively
  double nowpi, nowll; // log prior, log likelihood
  VectorXd nowtriDeltaQD; //precomputed components of the Wishart log likelihood, one for each microsatellite locus
  VectorXi nowqColors; // mapping that indicates which q tiles each vertex/deme falls into
  VectorXi nowmColors; // mapping that indicates which m tiles each vertex/deme falls into
  VectorXd nowW; // the within demes component, one entry for each deme in the habitat
  MatrixXd nowB; // the between demes component, one entry for each pair of demes in the habitat

  // Fixed constants -- necessary to scale diploid and haploid species slightly differently
  double ll_atfixdf, Wconst, Bconst;

  // Variables to store the results in :
  // Fixed size:
  MatrixXd mcmcmhyper;
  MatrixXd mcmcqhyper;
  MatrixXd mcmcthetas;
  MatrixXd mcmcpilogl;
  VectorXd mcmcmtiles;
  VectorXd mcmcqtiles;  
  // Variable length:
  vector<double> mcmcmRates;
  vector<double> mcmcqRates;
  vector<double> mcmcxCoord;
  vector<double> mcmcyCoord;
  vector<double> mcmcwCoord;
  vector<double> mcmczCoord;
  
  double EEMS_wishpdfln(const MatrixXd &B, const VectorXd &q, const VectorXd &sigma2, VectorXd &triDeltaQD) const;
  
};

#endif
