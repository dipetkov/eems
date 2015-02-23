#pragma once

#include "util.hpp"


class MCMC {
public:

  MCMC(const Params &params);
  ~MCMC();

  int currIter;
  int currStep;
  int numMCMCIter;
  int numBurnIter;
  int numThinIter;
  bool iterDone;
  bool finished;

  void start_iteration( );
  void end_iteration( );
  void add_to_okay_moves(const int type);
  void add_to_total_moves(const int type);
  int num_iters_to_save( ) const;
  int to_save_iteration( ) const;
  void change_update( );
  
  friend ostream& operator<<(ostream& out, const MCMC& mcmc);
    
private:

  int numTypes;
  vector<double> okayMoves;
  vector<double> totalMoves;

};
