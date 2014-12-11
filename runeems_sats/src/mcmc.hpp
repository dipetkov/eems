#pragma once

#include "util.hpp"


class MCMC {
public:

  MCMC(const Params &params);
  ~MCMC();

  int qTile;
  int mTile;
  int currIter;
  int currType;
  int numMCMCIter;
  int numBurnIter;
  int numThinIter;
  bool iterDone;
  bool finished;

  void start_iteration( );
  void end_iteration( );
  void addToOkayMoves( );
  void addToTotalMoves( );
  int num_iters_to_save( ) const;
  int to_save_iteration( ) const;

  void change_update(const int qtiles, const int mtiles);
  void output_proportions(ostream &out) const;
    
private:

  int numTypes;
  vector<double> okayMoves;
  vector<double> totalMoves;

};
