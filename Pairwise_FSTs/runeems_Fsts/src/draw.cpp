
#include "draw.hpp"

Draw::Draw( ) { }
Draw::~Draw( ) { }
void Draw::initialize(const long seed) {
  this->seed = seed;
  randgen = boost::mt19937(seed);
}
double Draw::runif( ) {
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
    runif(randgen, boost::uniform_real<>(0.0,1.0));
  return (randraw(runif));
}
int Draw::runif_int(const int min, const int max) {
  boost::variate_generator<boost::mt19937&, boost::uniform_int<> >
    runif(randgen, boost::uniform_int<>(min,max));
  return (randraw(runif));
}
double Draw::rnorm(const double mu, const double var) {
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >
    rnorm(randgen, boost::normal_distribution<>(mu,sqrt(var)));
  return (randraw(rnorm));
}
// Since rates are parameterized on the log scale, the support of
// the truncated normal distribution is symmetric about 0 in all cases,
// so instead of [xmin,xmax] it is sufficient to use [-bnd,+bnd]
double Draw::rtrnorm(const double mu, const double var, const double bnd) {
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >
    rnorm(randgen, boost::normal_distribution<>(mu,sqrt(var)));
  double x = -bnd - 1.0;
  while ((x< -bnd) || (x>bnd)) { x = randraw(rnorm); }
  return (x);
}  
double Draw::rinvgam(const double shape, const double scale) {
  boost::variate_generator<boost::mt19937&, boost::gamma_distribution<> >
    rgamma(randgen, boost::gamma_distribution<>(shape,1.0/scale));
  return (1.0/randraw(rgamma));
}
int Draw::rnegbin(const int r, const double p) {
  boost::variate_generator<boost::mt19937&, boost::random::negative_binomial_distribution<> >
    rnegbin(randgen, boost::random::negative_binomial_distribution<>(r,1.0-p));
  int k = 0;
  while (!k) { k = (int)randraw(rnegbin); }
  return (k);
}
