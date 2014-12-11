
#include "draw.hpp"

Draw::Draw( ) { }
Draw::~Draw( ) { }
void Draw::initialize(const long seed) {
  this->seed = seed;
  randgen = boost::mt19937(seed);
}
long Draw::get_seed( ) const {
  return (seed);
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
double Draw::rtnorm(const double mu, const double var, const double lob, const double upb) {
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >
    rnorm(randgen, boost::normal_distribution<>(mu,sqrt(var)));
  double x = lob - 1.0;
  while ((x<lob)||(x>upb)) { x = randraw(rnorm); }
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
