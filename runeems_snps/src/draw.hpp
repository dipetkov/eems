#pragma once

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/negative_binomial_distribution.hpp>

#include "util.hpp"

#ifndef DRAW_H
#define DRAW_H

/* Random number generation with Boost */

class Draw {
public:

  Draw( );
  ~Draw( );
  long seed;
  
  void initialize(const long seed);
  double runif( );
  int rnegbin(const int r, const double p);
  int runif_int(const int min, const int max);
  double rnorm(const double mu, const double var);
  double rtrnorm(const double mu, const double var, const double bnd);
  double rinvgam(const double shape, const double scale);

private:

  boost::mt19937 randgen;
  template<class T>
  double randraw(T &generator) { return generator( ); }
  
};

#endif
