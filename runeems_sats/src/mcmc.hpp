#pragma once

#include "util.hpp"

#ifndef MCMC_H
#define MCMC_H

enum MoveType {
  Q_VORONOI_RATE_UPDATE,
  Q_VORONOI_POINT_MOVE,
  Q_VORONOI_BIRTH_DEATH,
  M_VORONOI_RATE_UPDATE,
  M_MEAN_RATE_UPDATE,
  M_VORONOI_POINT_MOVE,
  M_VORONOI_BIRTH_DEATH,
  DF_UPDATE,
  UNKNOWN_MOVE_TYPE
};

class MCMC {
public:

  MCMC(const Params &params);
  ~MCMC();

  int currIter;
  int numMCMCIter;
  int numBurnIter;
  int numThinIter;
  bool finished;

  void end_iteration( );
  void add_to_okay_moves(const int type);
  void add_to_total_moves(const int type);
  int num_iters_to_save( ) const;
  int to_save_iteration( ) const;
  
  friend ostream& operator<<(ostream& out, const MCMC& mcmc);
    
private:

  int numTypes;
  vector<double> okayMoves;
  vector<double> totalMoves;

};

#endif
