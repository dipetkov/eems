#pragma once

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/negative_binomial_distribution.hpp>

#include "util.hpp"

/* Random number generation with Boost */

class Draw {
public:

  Draw( );
  ~Draw( );

  long seed;
  long get_seed( ) const;
  
  void initialize(const long seed);
  double runif( );
  int runif_int(const int min, const int max);
  double rnorm(const double mu, const double var);
  double rtnorm(const double mu, const double var, const double lob, const double upb);
  double rinvgam(const double shape, const double scale);
  int rnegbin(const int r, const double p);

private:

  boost::mt19937 randgen;
  template<class T>
  double randraw(T &generator) { return generator( ); }
  
};
